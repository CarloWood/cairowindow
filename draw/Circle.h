#pragma once

#include "Shape.h"
#include "cairowindow/LayerRegion.h"
#include "cairowindow/Line.h"
#include "cairowindow/Color.h"

namespace cairowindow::draw {

// List the additional members of CircleStyle (none).
#define cairowindow_Circle_FOREACH_MEMBER(X, ...)

// CircleStyle is derived from ShapeStyle.
#define cairowindow_Circle_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Shape_FOREACH_STYLE_MEMBER(X, __VA_ARGS__) \
  cairowindow_Circle_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for CircleStyle.
struct CircleStyleParamsDefault : ShapeStyleParamsDefault
{
  // Override defaults from ShapeStyleParamsDefault.
  static constexpr Color line_color = color::black;
  static constexpr double line_width = 1.0;
  static constexpr Color fill_color = color::transparent;
  static constexpr ShapePosition position = at_corner;          // The center of the circle is normally at the top-left corner.
  static constexpr ShapeEnum shape = ellipse;
};

// Declare CircleStyle, derived from ShapeStyle.
DECLARE_STYLE_WITH_BASE(Circle, Shape, CircleStyleParamsDefault);

class Circle : public Shape
{
 public:
  Circle(cairowindow::Rectangle const& geometry, CircleStyle style) : Shape(geometry, style({.shape = ellipse})) { }

  CircleStyle const& style() const
  {
    // The static cast is safe because (currently) CircleStyle does not have additional members.
    return static_cast<CircleStyle const&>(style_);
  }
};

} // namespace cairowindow::draw
