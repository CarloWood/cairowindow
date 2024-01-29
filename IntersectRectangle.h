#pragma once

#include "Rectangle.h"
#include <algorithm>

namespace cairowindow {

class StrokeExtents;

class IntersectRectangle
{
 private:
  double x1_;
  double y1_;
  double x2_;
  double y2_;

 public:
  IntersectRectangle() = default;
  inline IntersectRectangle(StrokeExtents const& stroke_extents);

  IntersectRectangle(Rectangle const& rectangle) :
    x1_(rectangle.offset_x()), y1_(rectangle.offset_y()),
    x2_(rectangle.offset_x() + rectangle.width()), y2_(rectangle.offset_y() + rectangle.height()) { }

  double x1() const { return x1_; }
  double x2() const { return x2_; }
  double y1() const { return y1_; }
  double y2() const { return y2_; }

  IntersectRectangle(IntersectRectangle rect1, IntersectRectangle rect2) :
    x1_(std::max(rect1.x1_, rect2.x1_)), y1_(std::max(rect1.y1_, rect2.y1_)),
    x2_(std::min(rect1.x2_, rect2.x2_)), y2_(std::min(rect1.y2_, rect2.y2_)) { }
};

} // namespace cairowindow
