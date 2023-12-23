#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#include "cairowindow/Defaults.h"
#include <string>

namespace cairowindow::draw {

enum TextPosition
{
  above_right_off,
  above_left_off,
  below_right_off,
  below_left_off,
  centered_right_off,
  centered_left_off,
  centered_above,
  centered_below,
  centered
};

DECLARE_DEFAULTS_HAS_MEMBER(position)
DECLARE_DEFAULTS_HAS_MEMBER(font_size)
DECLARE_DEFAULTS_HAS_MEMBER(color)
DECLARE_DEFAULTS_HAS_MEMBER(font_family)
DECLARE_DEFAULTS_HAS_MEMBER(offset)
DECLARE_DEFAULTS_HAS_MEMBER(rotation)

struct DefaultTextStyleDefaults
{
  static constexpr draw::TextPosition position = draw::above_right_off;
  static constexpr double font_size = 12.0;
  static constexpr Color color = color::black;
  static constexpr char const* font_family = "sans-serif";
  static constexpr double offset = 10.0;
  static constexpr double rotation = 0.0;        // Clock-wise rotation in radians.
};

template<typename Defaults = DefaultTextStyleDefaults>
struct TextStyle
{
  TextPosition position   = DEFAULT_FROM(Defaults, position);
  double font_size        = DEFAULT_FROM(Defaults, font_size);
  Color color             = DEFAULT_FROM(Defaults, color);
  std::string font_family = DEFAULT_FROM(Defaults, font_family);
  double offset           = DEFAULT_FROM(Defaults, offset);
  double rotation         = DEFAULT_FROM(Defaults, rotation);

  void setup(cairo_t* cr)
  {
    cairo_select_font_face(cr, font_family.c_str(), CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, font_size);
    cairo_set_source_rgb(cr, color.red(), color.green(), color.blue());
  }
};

struct TextStyleNoDefault : TextStyle<>
{
  TextStyleNoDefault() = default;

  template<typename Defaults>
  TextStyleNoDefault(TextStyle<Defaults> const& text_style) :
    TextStyle<>{
      .position = text_style.position,
      .font_size = text_style.font_size,
      .color = text_style.color,
      .font_family = text_style.font_family,
      .offset = text_style.offset,
      .rotation = text_style.rotation}
  { }
};

class Text : public LayerRegion
{
 private:
  TextStyleNoDefault style_;
  std::string text_;
  double pos_x_{};
  double pos_y_{};

 public:
  Text() = default;
  Text(std::string const& text, double pos_x, double pos_y, TextStyleNoDefault style) :
    style_(style), text_(text), pos_x_(pos_x), pos_y_(pos_y) { }

  bool is_defined() const
  {
    return !text_.empty();
  }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::notice, "Text::do_draw(cr)");
    style_.setup(cr);
    cairo_text_extents_t extents;
    cairo_text_extents(cr, text_.c_str(), &extents);
    cairo_translate(cr, pos_x_, pos_y_);
    cairo_rotate(cr, style_.rotation);
    double tx = 0;
    double ty = 0;
    switch (style_.position)
    {
      case above_right_off:
        break;
      case above_left_off:
        tx = -extents.x_advance;
        break;
      case below_right_off:
        ty = -extents.y_bearing;
        break;
      case below_left_off:
        tx = -extents.x_advance;
        ty = -extents.y_bearing;
        break;
      case centered_right_off:
        ty = -0.5 * extents.y_bearing;
        break;
      case centered_left_off:
        tx = -extents.x_advance;
        ty = -0.5 * extents.y_bearing;
        break;
      case centered_above:
        tx = -extents.x_bearing - 0.5 * extents.width;
        break;
      case centered_below:
        tx = -extents.x_bearing - 0.5 * extents.width;
        ty = -extents.y_bearing;
        break;
      case centered:
        tx = -extents.x_bearing - 0.5 * extents.width;
        ty = -0.5 * extents.y_bearing;
    }
    cairo_translate(cr, tx, ty);
    cairo_show_text(cr, text_.c_str());
    cairo_stroke(cr);
    double x1 = extents.x_bearing;
    double y1 = extents.y_bearing;
    double x2 = extents.x_bearing + extents.width;
    double y2 = extents.y_bearing + extents.height;
    cairo_matrix_t m;
    cairo_get_matrix(cr, &m);
    cairo_matrix_transform_point(&m, &x1, &y1);
    cairo_matrix_transform_point(&m, &x2, &y2);
    return { x1, y1, x2, y2 };
  }
};

} // namespace cairowindow::draw
