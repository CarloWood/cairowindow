#pragma once

#include "cairowindow/Point.h"
#include "cairowindow/Direction.h"
#include "Shape.h"
#include <array>

namespace cairowindow::draw {

// List the additional members of ArrowHeadBaseStyle (none).
#define cairowindow_ArrowHeadBase_FOREACH_MEMBER(X, ...)

// ArrowHeadBaseStyle is derived from ShapeStyle.
#define cairowindow_ArrowHeadBase_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Shape_FOREACH_STYLE_MEMBER(X, __VA_ARGS__) \
  cairowindow_ArrowHeadBase_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for ArrowHeadBaseStyle.
struct ArrowHeadBaseStyleParamsDefault : ShapeStyleParamsDefault
{
};

// Declare ArrowHeadBaseStyle, derived from ShapeStyle.
DECLARE_STYLE_WITH_BASE(ArrowHeadBase, Shape, ArrowHeadBaseStyleParamsDefault);

// Extend ArrowHeadStyle with a member function.
class ArrowHeadStyle : public ArrowHeadBaseStyle
{
 public:
  using ArrowHeadBaseStyle::ArrowHeadBaseStyle;
  ArrowHeadStyle(ArrowHeadBaseStyle const& style) : ArrowHeadBaseStyle(style) { }

  // Return an index usable for s_arrow_head_size.
  int arrow() const { return m_shape - none_arrow_shape; }
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
        s_arrow_head_size[style.arrow()].width, s_arrow_head_size[style.arrow()].height}, style,
        direction.as_angle()),
    tip_x_(tip_x), tip_y_(tip_y), direction_(direction), style_(style) { }
};

} // namespace cairowindow::draw
