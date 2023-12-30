#pragma once

#include "Point.h"
#include "utils/square.h"
#include <cmath>

namespace cairowindow {

class LinePiece
{
 private:
  Point from_;
  Point to_;

 public:
  LinePiece(Point const& from, Point const& to) : from_(from), to_(to) { }

  Point const& from() const { return from_; }
  Point const& to() const { return to_; }
  double length() const { return std::sqrt(utils::square(from_.x() - to_.x()) + utils::square(from_.y() - to_.y())); };
};

} // namespace cairowindow
