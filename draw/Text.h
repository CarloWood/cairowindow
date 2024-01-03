#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#include "cairowindow/Defaults.h"
#include <string>
#ifdef CWDEBUG
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

enum TextPosition
{
  undefined,
  above_right_of,
  above_left_of,
  below_right_of,
  below_left_of,
  centered_right_of,
  centered_left_of,
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
  static constexpr draw::TextPosition position = draw::above_right_of;
  static constexpr double font_size = 12.0;
  static constexpr Color color = color::black;
  static constexpr char const* font_family = "sans-serif";
  static constexpr double offset = 0.0;
  static constexpr double rotation = 0.0;        // Clock-wise rotation in radians.
};

struct TextStyleDelta
{
  static constexpr double undefined_offset_magic = 12345678.9;
  static constexpr double undefined_rotation_magic = -3000000.0;

  TextPosition position   = undefined;
  double font_size        = -1.0;
  Color color{};
  std::string font_family{};
  double offset           = undefined_offset_magic;
  double rotation         = undefined_rotation_magic;
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
    DoutEntering(dc::notice, "TextStyle<" << libcwd::type_info_of<Defaults>().demangled_name() << ">::setup(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_select_font_face(cr, font_family.c_str(), CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, font_size);
    cairo_set_source_rgb(cr, color.red(), color.green(), color.blue());
  }

  TextStyle operator()(TextStyleDelta delta)
  {
    TextStyle result{*this};
    if (delta.position != undefined)
      result.position = delta.position;
    if (delta.font_size != -1.0)
      result.font_size = delta.font_size;
    if (delta.color.is_defined())
      result.color = delta.color;
    if (!delta.font_family.empty())
      result.font_family = delta.font_family;
    if (delta.offset != TextStyleDelta::undefined_offset_magic)
      result.offset = delta.offset;
    if (delta.rotation != TextStyleDelta::undefined_rotation_magic)
      result.rotation = delta.rotation;
    return result;
  }

#ifdef CWDEBUG
  friend std::ostream& operator<<(std::ostream& os, TextStyle const* text_style_ptr)
  {
    os << "TextStyle*";
    return os;
  }
#endif
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

  // Accessors.
  TextStyleNoDefault const& style() const { return style_; }
  std::string const& text() const { return text_; }
  double pos_x() const { return pos_x_; }
  double pos_y() const { return pos_y_; }

  // Allow moving the text after construction (but before adding it to a layer!).
  void rel_move_to(double delta_x, double delta_y)
  {
    pos_x_ += delta_x;
    pos_y_ += delta_y;
  }

  void move_to(double x, double y)
  {
    pos_x_ = x;
    pos_y_ = y;
  }

  bool is_defined() const
  {
    return !text_.empty();
  }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::notice, "Text::do_draw(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    style_.setup(cr);
    cairo_text_extents_t extents;
    Dout(dc::notice, "Drawing \"" << text_ << "\" at position (" << pos_x_ << ", " << pos_y_ << ")");
    cairo_text_extents(cr, text_.c_str(), &extents);
    cairo_translate(cr, pos_x_, pos_y_);
    cairo_rotate(cr, style_.rotation);
    double tx = 0;
    double ty = 0;
    switch (style_.position)
    {
      case undefined:
        ASSERT(false);
      case above_right_of:
        tx += style_.offset;
        break;
      case above_left_of:
        tx -= style_.offset;
        tx -= extents.x_advance;
        break;
      case below_right_of:
        tx += style_.offset;
        ty -= extents.y_bearing;
        break;
      case below_left_of:
        tx -= style_.offset;
        tx -= extents.x_advance;
        ty -= extents.y_bearing;
        break;
      case centered_right_of:
        tx += style_.offset;
        ty -= 0.5 * extents.y_bearing;
        break;
      case centered_left_of:
        tx -= style_.offset;
        tx -= extents.x_advance;
        ty -= 0.5 * extents.y_bearing;
        break;
      case centered_above:
        ty -= style_.offset;
        tx -= extents.x_bearing + 0.5 * extents.width;
        break;
      case centered_below:
        ty += style_.offset;
        tx -= extents.x_bearing + 0.5 * extents.width;
        ty -= extents.y_bearing;
        break;
      case centered:
        tx -= extents.x_bearing + 0.5 * extents.width;
        ty -= 0.5 * extents.y_bearing;
        break;
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
    // Due to (90 degree) rotation the coordinates are possibly no longer ordered.
    // The second coordinate must be the largest.
    if (x2 < x1)
      std::swap(x1, x2);
    if (y2 < y1)
      std::swap(y1, y2);
    Dout(dc::notice, "Returning (" << x1 << ", " << y1 << ", " << x2 << ", " << y2 << ")");
    return { x1, y1, x2, y2 };
  }
};

} // namespace cairowindow::draw
