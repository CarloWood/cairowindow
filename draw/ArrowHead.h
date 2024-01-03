#pragma once

#include "cairowindow/Point.h"
#include "cairowindow/Direction.h"
#include "Shape.h"
#include <array>

namespace cairowindow::draw {

struct ArrowHeadStyle
{
  // Must start with the same style variables as ShapeStyle, because we do a reinterpret_cast between the two!
  Color line_color = color::transparent;
  Color fill_color = color::transparent;
  double line_width = 1.0;
  ShapePosition position = at_tip;
  ShapeEnum shape = open_arrow_shape;

  // Return an index usable for s_arrow_head_size.
  int arrow() const { return shape - none_arrow_shape; }
};

class ArrowHead : public Shape
{
 private:
  double tip_x_;
  double tip_y_;
  Direction direction_;
  ArrowHeadStyle style_;

  struct Size { double width; double height; };         // In pixels.
  static std::array<Size, number_of_arrow_shapes> s_arrow_head_size;

 public:
  ArrowHead(double tip_x, double tip_y, Direction direction, ArrowHeadStyle style) :
    Shape({tip_x - s_arrow_head_size[style.arrow()].width, tip_y - 0.5 * s_arrow_head_size[style.arrow()].height,
        s_arrow_head_size[style.arrow()].width, s_arrow_head_size[style.arrow()].height}, reinterpret_cast<ShapeStyle const&>(style),
        direction.as_angle()),
    tip_x_(tip_x), tip_y_(tip_y), direction_(direction), style_(style) { }
};

} // namespace cairowindow::draw
