#pragma once

#include "Point.h"

namespace cairowindow::cs {

template<CS cs>
class Circle
{
 private:
  Point<cs> center_;
  double radius_{};

 public:
  Circle() = default;
  Circle(Point<cs> const& center, double radius) : center_(center), radius_(radius) { }

  Point<cs> const& center() const { return center_; }
  double radius() const { return radius_; }
};

} // namespace cairowindow::cs

