#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"

namespace cairowindow::draw {

enum ShapeEnum
{
  rectangle,
  ellipse,
  triangle,
  triangle_up = triangle,
  triangle_down,
  triangle_left,
  triangle_right,
  rhombus,
  diamond,
  cross,
  plus,
  star,
  // Shape is used to also cover arrow heads.
  none_arrow_shape,
  open_arrow_shape,
  filled_arrow_shape,
  diamond_arrow_shape,
  circle_arrow_shape
};

static constexpr int number_of_shapes = star + 1;
static constexpr int number_of_arrow_shapes = circle_arrow_shape + 1 - none_arrow_shape;

ShapeEnum next_shape();

enum ShapePosition
{
  at_center,
  at_corner,
  at_tip
};

struct ShapeStyle
{
  // Must start with the same style variables as ArrowHeadStyle and CircleStyle, because we do a reinterpret_cast between these!
  Color line_color = color::transparent;
  Color fill_color = color::transparent;
  double line_width = 2.0;
  ShapePosition position = at_center;
  // This only exists in ShapeStyle and ArrowHeadStyle.
  ShapeEnum shape = rectangle;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

class Shape : public LayerRegion
{
 protected:
  Rectangle geometry_;
  ShapeStyle style_;
  double rotation_;

 public:
  Shape(Rectangle const& geometry, ShapeStyle style, double rotation = {}) : geometry_(geometry), style_(style), rotation_(rotation) { }

  ShapeStyle& style() { return style_; }
  ShapeStyle const& style() const { return style_; }

 private:
  StrokeExtents do_draw(cairo_t* cr) override;
};

} // namespace cairowindow::draw
