#pragma once

#include "math/cs/Point.h"

namespace cairowindow::cs {

template<CS cs>
class Circle
{
 private:
  math::cs::Point<cs> center_;
  double radius_{};

 public:
  Circle() = default;
  Circle(math::cs::Point<cs> const& center, double radius) : center_(center), radius_(radius) { }

  math::cs::Point<cs> const& center() const { return center_; }
  double radius() const { return radius_; }
};

} // namespace cairowindow::cs

