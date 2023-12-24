#pragma once

#include "Shape.h"
#include <array>

namespace cairowindow::draw {

class Point
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

 private:
  Shape shape_;

 public:
  Point(double x, double y, int color_index, int filled_shape) :
    shape_({x, y, is_filled(filled_shape) ? 5.0 : 4.0, is_filled(filled_shape) ? 5.0 : 4.0}, {
        .line_color = is_filled(filled_shape) ? color::transparent : Color::get_color(color_index),
        .fill_color = is_filled(filled_shape) ? Color::get_color(color_index) : color::transparent,
        .shape = get_shape(filled_shape), .at_corner = true })
  {
  }

  operator LayerRegion*()
  {
    return &shape_;
  }

 private:
  bool is_filled(int filled_shape)
  {
    return filled_shapes[filled_shape % number_of_shapes].filled;
  }

  ShapeEnum get_shape(int filled_shape)
  {
    return filled_shapes[filled_shape % number_of_shapes].shape;
  }
};

} // namespace cairowindow::draw

