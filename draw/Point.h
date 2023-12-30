#pragma once

#include "Shape.h"
#include <array>

namespace cairowindow::draw {

struct PointStyle {
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

 private:
  int color_index_;
  int filled_shape_;

 public:
  PointStyle(int color_index, int filled_shape) : color_index_(color_index), filled_shape_(filled_shape) { }

  bool is_filled() const
  {
    return filled_shapes[filled_shape_ % number_of_shapes].filled;
  }

  ShapeEnum get_shape() const
  {
    return filled_shapes[filled_shape_ % number_of_shapes].shape;
  }

  Color line_color() const
  {
    return is_filled() ? color::transparent : Color::get_color(color_index_);
  };

  Color fill_color() const
  {
    return is_filled() ? Color::get_color(color_index_) : color::transparent;
  }
};

class Point : public Shape
{
 public:
  Point(double x, double y, PointStyle style) :
    Shape({x, y, style.is_filled() ? 5.0 : 4.0, style.is_filled() ? 5.0 : 4.0},
          { .line_color = style.line_color(), .fill_color = style.fill_color(), .at_corner = true, .shape = style.get_shape() })
  {
    DoutEntering(dc::notice, "Point(" << x << ", " << y << ", style) [" << this << "]");
  }

  ~Point()
  {
    DoutEntering(dc::notice, "~Point() [" << this << "]");
  }
};

} // namespace cairowindow::draw

