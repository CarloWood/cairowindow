#pragma once

#include "Direction.h"
#include "Point.h"

namespace cairowindow {

class Line
{
 private:
  Direction normal_;
  Point point_;

 public:
  Line(Direction const& normal, Point const& point) : normal_(normal), point_(point) { }

  Direction normal() const { return normal_; }
  Point point() const { return point_; }
};

} // namespace cairowindow
