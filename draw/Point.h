#pragma once

#include "Shape.h"
#include <array>
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#endif

namespace cairowindow::draw {

struct PointStyleDelta
{
  static constexpr int undefined_magic = -1;

  int color_index = undefined_magic;
  int filled_shape = undefined_magic;
};

struct PointStyle
{
 public:
  struct FilledShape
  {
    ShapeEnum shape;
    bool filled;
  };

  static constexpr int number_of_shapes = 15;
  static constexpr std::array<FilledShape, number_of_shapes> filled_shapes = {{
    { rectangle, false },
    { plus, false },
    { triangle, true },
    { ellipse, true },
    { rhombus, false },
    { cross, false },
    { diamond, true, },
    { triangle_down, false },
    { rectangle, true },
    { triangle, false },
    { ellipse, false },
    { rhombus, true },
    { diamond, false },
    { triangle_down, true },
    { star, false }
  }};

 public:
  int color_index;
  int filled_shape;

  bool is_filled() const
  {
    return filled_shapes[filled_shape % number_of_shapes].filled;
  }

  ShapeEnum get_shape() const
  {
    return filled_shapes[filled_shape % number_of_shapes].shape;
  }

  Color line_color() const
  {
    return is_filled() ? color::transparent : Color::get_color(color_index);
  };

  Color fill_color() const
  {
    return is_filled() ? Color::get_color(color_index) : color::transparent;
  }

  PointStyle operator()(PointStyleDelta delta)
  {
    PointStyle result{*this};
    if (delta.color_index != PointStyleDelta::undefined_magic)
      result.color_index = delta.color_index;
    if (delta.filled_shape != PointStyleDelta::undefined_magic)
      result.filled_shape = delta.filled_shape;
    return result;
  }
};

class Point : public Shape
{
 private:
  PointStyle point_style_;

 public:
  Point(double x, double y, PointStyle style) :
    Shape({x, y, style.is_filled() ? 5.0 : 4.0, style.is_filled() ? 5.0 : 4.0},
          { .line_color = style.line_color(), .fill_color = style.fill_color(), .position = at_corner, .shape = style.get_shape() }),
    point_style_(style)
  {
    DoutEntering(dc::cairowindow, "Point(" << x << ", " << y << ", style) [" << this << "]");
  }

  ~Point()
  {
    DoutEntering(dc::cairowindow, "~Point() [" << this << "]");
  }

  PointStyle const& point_style() const { return point_style_; }
};

} // namespace cairowindow::draw

