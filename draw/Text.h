#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#include "cairowindow/Style.h"
#include <string>
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

enum TextPosition
{
  undefined_text_position,
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

// List the members of TextStyle.
#define cairowindow_TextBase_FOREACH_MEMBER(X, ...) \
  X(draw::TextPosition, position, undefined_text_position, __VA_ARGS__) \
  X(double, font_size, -1.0, __VA_ARGS__) \
  X(Color, color, Color{}, __VA_ARGS__) \
  X(std::string, font_family, "", __VA_ARGS__) \
  X(double, offset, 12345678.9, __VA_ARGS__) \
  X(double, rotation, -3000000.0, __VA_ARGS__)

// Mandatory macro.
#define cairowindow_TextBase_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_TextBase_FOREACH_MEMBER(X, __VA_ARGS__)

struct TextStyleParamsDefault
{
  static constexpr draw::TextPosition position = draw::above_right_of;
  static constexpr double font_size = 12.0;
  static constexpr Color color = color::black;
  static constexpr char const* font_family = "sans-serif";
  static constexpr double offset = 0.0;
  static constexpr double rotation = 0.0;        // Clock-wise rotation in radians.
};

// Declare TextBaseStyle.
DECLARE_STYLE(TextBase, TextStyleParamsDefault);

// Extend TextStyleBase with a member function.
class TextStyle : public TextBaseStyle
{
 public:
  using TextBaseStyle::TextBaseStyle;
  TextStyle(TextBaseStyle const& text_style) : TextBaseStyle(text_style) { }

  void setup(cairo_t* cr)
  {
    DoutEntering(dc::cairowindow, "TextStyle::setup(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_select_font_face(cr, m_font_family.c_str(), CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, m_font_size);
    cairo_set_source_rgb(cr, m_color.red(), m_color.green(), m_color.blue());
  }
};

class Text : public LayerRegion
{
 private:
  TextStyle style_;
  std::string text_;
  double pos_x_{};
  double pos_y_{};

 public:
  Text(std::string const& text, double pos_x, double pos_y, TextStyle style) :
    style_(style), text_(text), pos_x_(pos_x), pos_y_(pos_y) { }

  // Accessors.
  TextStyle const& style() const { return style_; }
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
    DoutEntering(dc::cairowindow, "Text::do_draw(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    style_.setup(cr);
    cairo_text_extents_t extents;
    Dout(dc::cairowindow, "Drawing \"" << text_ << "\" at position (" << pos_x_ << ", " << pos_y_ << ")");
    cairo_text_extents(cr, text_.c_str(), &extents);
    cairo_translate(cr, pos_x_, pos_y_);
    cairo_rotate(cr, style_.rotation());
    double tx = 0;
    double ty = 0;
    switch (style_.position())
    {
      case undefined_text_position:
        ASSERT(false);
      case above_right_of:
        tx += style_.offset();
        break;
      case above_left_of:
        tx -= style_.offset();
        tx -= extents.x_advance;
        break;
      case below_right_of:
        tx += style_.offset();
        ty -= extents.y_bearing;
        break;
      case below_left_of:
        tx -= style_.offset();
        tx -= extents.x_advance;
        ty -= extents.y_bearing;
        break;
      case centered_right_of:
        tx += style_.offset();
        ty -= 0.5 * extents.y_bearing;
        break;
      case centered_left_of:
        tx -= style_.offset();
        tx -= extents.x_advance;
        ty -= 0.5 * extents.y_bearing;
        break;
      case centered_above:
        ty -= style_.offset();
        tx -= extents.x_bearing + 0.5 * extents.width;
        break;
      case centered_below:
        ty += style_.offset();
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
    Dout(dc::cairowindow, "Returning (" << x1 << ", " << y1 << ", " << x2 << ", " << y2 << ")");
    return { x1, y1, x2, y2 };
  }
};

} // namespace cairowindow::draw
