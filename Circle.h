#pragma once

#include "Point.h"

namespace cairowindow {

class Circle
{
 private:
  Point center_;
  double radius_;

 public:
  Circle() = default;
  Circle(Point const& center, double radius) : center_(center), radius_(radius) { }

  Point const& center() const { return center_; }
  double radius() const { return radius_; }
};

} // namespace cairowindow

