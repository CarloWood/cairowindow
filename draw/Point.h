#pragma once

#include "Shape.h"
#include "cairowindow/Style.h"
#include <array>
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#endif

namespace cairowindow::draw {

#define cairowindow_PointBase_FOREACH_MEMBER(X, ...) \
  X(int, color_index, -1, __VA_ARGS__) \
  X(int, filled_shape, -1, __VA_ARGS__)

#define cairowindow_PointBase_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_PointBase_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for PointStyle.
struct PointStyleParamsDefault
{
  static constexpr int color_index = 0;
  static constexpr int filled_shape = 0;
};

// Declare PointStyle.
DECLARE_STYLE(PointBase, PointStyleParamsDefault);

// Extend PointBaseStyle with some member functions.
struct PointStyle : public PointBaseStyle
{
 public:
  using PointBaseStyle::PointBaseStyle;
  PointStyle(PointBaseStyle const& point_style) : PointBaseStyle(point_style) { }

  struct FilledShape
  {
    ShapeEnum shape;
    bool filled;
  };

  static constexpr int number_of_shapes = 16;
  static constexpr std::array<FilledShape, number_of_shapes> filled_shapes = {{
    { rectangle, false },               // 0
    { plus, false },
    { triangle, true },
    { ellipse, true },
    { rhombus, false },
    { cross, false },                   // 5
    { diamond, true, },
    { triangle_down, false },
    { rectangle, true },
    { triangle, false },
    { ellipse, false },                 // 10
    { rhombus, true },
    { diamond, false },
    { triangle_down, true },
    { star, false },
    { dot, false }                      // 15
  }};

 public:
  bool is_filled() const
  {
    return filled_shapes[m_filled_shape % number_of_shapes].filled;
  }

  ShapeEnum get_shape() const
  {
    return filled_shapes[m_filled_shape % number_of_shapes].shape;
  }

  Color line_color() const
  {
    return is_filled() ? color::transparent : Color::get_color(m_color_index);
  };

  Color fill_color() const
  {
    return is_filled() ? Color::get_color(m_color_index) : color::transparent;
  }
};

class Point : public Shape
{
 private:
  PointStyle point_style_;

 public:
  Point(double x, double y, PointStyle style) :
    Shape({x, y, style.is_filled() ? 5.0 : 4.0, style.is_filled() ? 5.0 : 4.0},
          ShapeStyleParams{ .line_color = style.line_color(), .fill_color = style.fill_color(), .position = at_corner, .shape = style.get_shape() }),
    point_style_(style)
  {
    DoutEntering(dc::cairowindow, "Point(" << x << ", " << y << ", style) [" << this << "]");
  }

  Point(cairowindow::cs::Point<CS::pixels> const& point, PointStyle style) : Point(point.x(), point.y(), style) { }

  ~Point()
  {
    DoutEntering(dc::cairowindow, "~Point() [" << this << "]");
  }

  PointStyle const& point_style() const { return point_style_; }
};

} // namespace cairowindow::draw

