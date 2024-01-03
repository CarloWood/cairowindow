#pragma once

#include "Direction.h"
#include "Point.h"

namespace cairowindow {

class Line
{
 private:
  Direction direction_;
  Point point_;

 public:
  Line(Direction const& direction, Point const& point) : direction_(direction), point_(point) { }

  Direction const& direction() const { return direction_; }
  Point const& point() const { return point_; }

  operator Direction const&() const { return direction_; }

  Point intersection_with(Line const& line2) const;
};

} // namespace cairowindow
