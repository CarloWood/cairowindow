#pragma once

#include "Rectangle.h"
#include "utils/has_print_on.h"
#include <algorithm>

namespace cairowindow {
using utils::has_print_on::operator<<;

class StrokeExtents;

template<CS cs>
requires (cs == CS::pixels || cs == CS::plot)
class IntersectRectangle
{
 private:
  double x1_;
  double y1_;
  double x2_;
  double y2_;

 public:
  IntersectRectangle() = default;

  inline IntersectRectangle(StrokeExtents const& stroke_extents) requires (cs == CS::pixels);   // Defined after StrokeExtents.

  IntersectRectangle(cs::Rectangle<cs> const& rectangle) :
    x1_(rectangle.offset_x()), y1_(rectangle.offset_y()),
    x2_(rectangle.offset_x() + rectangle.width()), y2_(rectangle.offset_y() + rectangle.height()) { }

  double x1() const { return x1_; }
  double x2() const { return x2_; }
  double y1() const { return y1_; }
  double y2() const { return y2_; }

  IntersectRectangle(IntersectRectangle rect1, IntersectRectangle rect2) :
    x1_(std::max(rect1.x1_, rect2.x1_)), y1_(std::max(rect1.y1_, rect2.y1_)),
    x2_(std::min(rect1.x2_, rect2.x2_)), y2_(std::min(rect1.y2_, rect2.y2_)) { }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{x1:" << x1_ << ", y1:" << y1_ << ", x2:" << x2_ << ", y2:" << y2_ << "}";
  }
#endif
};

} // namespace cairowindow
