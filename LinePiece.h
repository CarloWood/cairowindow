#pragma once

#include "Point.h"
#include "utils/square.h"
#include <cmath>
#include <memory>

namespace cairowindow {

class LinePiece
{
 private:
  Point from_;
  Point to_;

 public:
  LinePiece() = default;
  LinePiece(Point const& from, Point const& to) : from_(from), to_(to) { }

  Point const& from() const { return from_; }
  Point const& to() const { return to_; }
  double length() const { return std::sqrt(utils::square(from_.x() - to_.x()) + utils::square(from_.y() - to_.y())); };
};

namespace draw {
class Line;
} // namespace draw

namespace plot {
class Plot;

//--------------------------------------------------------------------------
// LinePiece

class LinePiece : public cairowindow::LinePiece
{
 public:
  using cairowindow::LinePiece::LinePiece;
  LinePiece(cairowindow::LinePiece const& line_piece) : cairowindow::LinePiece(line_piece) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

} // namespace plot
} // namespace cairowindow
