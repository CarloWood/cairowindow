#pragma once

#include "LineStyle.h"
#include "ShapePosition.h"
#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/Rectangle.h"
#include "utils/square.h"

namespace cairowindow::draw {

enum ShapeEnum
{
  undefined_shape,
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
  circle_arrow_shape,
  // Shape used for a single pixel.
  dot
};

static constexpr int number_of_shapes = star + 1;
static constexpr int number_of_arrow_shapes = circle_arrow_shape + 1 - none_arrow_shape;

ShapeEnum next_shape();

// List the additional members of ShapeStyle.
#define cairowindow_Shape_FOREACH_MEMBER(X, ...) \
  X(Color, fill_color, Color{}, __VA_ARGS__) \
  X(ShapePosition, position, undefined_shape_position, __VA_ARGS__) \
  X(ShapeEnum, shape, undefined_shape, __VA_ARGS__)

// ShapeStyle is derived from LineStyle.
#define cairowindow_Shape_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Line_FOREACH_STYLE_MEMBER(X, __VA_ARGS__) \
  cairowindow_Shape_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for ShapeStyle.
struct ShapeStyleParamsDefault : public LineStyleParamsDefault
{
  // Override defaults from CircleStyleParamsDefault.
  static constexpr Color line_color = color::transparent;
  // New defaults.
  static constexpr Color fill_color = color::transparent;
  static constexpr ShapePosition position = at_center;
  static constexpr ShapeEnum shape = rectangle;
};

// Declare ShapeStyle, derived from LineStyle.
DECLARE_STYLE_WITH_BASE(Shape, Line, ShapeStyleParamsDefault);

class Shape : public LayerRegion
{
 protected:
  cairowindow::Rectangle geometry_;
  ShapeStyle style_;
  double rotation_;
  double arrow_overshoot_;

 public:
  Shape(cairowindow::Rectangle const& geometry, ShapeStyle style, double rotation = {}) :
    geometry_(geometry), style_(style), rotation_(rotation)
    {
      switch (style.shape())
      {
        case open_arrow_shape:
        case filled_arrow_shape:
          arrow_overshoot_ = 0.5 * style_.line_width() * std::sqrt(::utils::square(2.0 * geometry_.width() / geometry_.height()) + 1.0);
          break;
        case diamond_arrow_shape:
          arrow_overshoot_ = 0.25 * style_.line_width() * std::sqrt(::utils::square(geometry_.width() / geometry_.height()) + 1.0);
          break;
        case circle_arrow_shape:
          arrow_overshoot_ = 0.5 * style_.line_width();
          break;
        default:        // Stop compiler from complaining.
          break;
      }
    }

  cairowindow::Rectangle const& geometry() const { return geometry_; }
  ShapeStyle const& style() const { return style_; }

  double arrow_overshoot() const { return arrow_overshoot_; }

 private:
  StrokeExtents do_draw(cairo_t* cr) override;
};

} // namespace cairowindow::draw
