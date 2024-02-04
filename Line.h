#pragma once

#include "Direction.h"
#include "Point.h"

namespace cairowindow {

class Line
{
 private:
  Point point_;
  Direction direction_;

 public:
  Line(Point const& point, Direction const& direction) : point_(point), direction_(direction) { }

  Point const& point() const { return point_; }
  Direction const& direction() const { return direction_; }

  operator Direction const&() const { return direction_; }

  Point intersection_with(Line const& line2) const;
};

} // namespace cairowindow
