#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "utils/square.h"

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
  double arrow_overshoot_;

 public:
  Shape(Rectangle const& geometry, ShapeStyle style, double rotation = {}) :
    geometry_(geometry), style_(style), rotation_(rotation)
    {
      switch (style.shape)
      {
        case open_arrow_shape:
        case filled_arrow_shape:
          arrow_overshoot_ = 0.5 * style_.line_width * std::sqrt(::utils::square(2.0 * geometry_.width() / geometry_.height()) + 1.0);
          break;
        case diamond_arrow_shape:
          arrow_overshoot_ = 0.25 * style_.line_width * std::sqrt(::utils::square(geometry_.width() / geometry_.height()) + 1.0);
          break;
        case circle_arrow_shape:
          arrow_overshoot_ = 0.5 * style_.line_width;
          break;
        default:        // Stop compiler from complaining.
          break;
      }
    }

  ShapeStyle& style() { return style_; }
  ShapeStyle const& style() const { return style_; }

  double arrow_overshoot() const { return arrow_overshoot_; }

 private:
  StrokeExtents do_draw(cairo_t* cr) override;
};

} // namespace cairowindow::draw
