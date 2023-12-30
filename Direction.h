#pragma once

#include "LinePiece.h"
#include <cmath>

namespace cairowindow {

class Line;

class Direction
{
 private:
  double x_;    // x,y is a unit vector pointing in the intended direction.
  double y_;

 public:
  // Construct a Direction that points in the direction theta (in radians): an angle with the positive x-axis.
  Direction(double theta) : x_(std::cos(theta)), y_(std::sin(theta)) { }

  // Construct a Direction pointing from origin to (x, y).
  Direction(double x, double y) : x_(x), y_(y) { }

  // Construct a Direction from two points. If the second point is not given it defaults to the origin.
  // The direction is from the second argument (or origin) to the first argument.
  Direction(Point const& to, Point const& from = Point{0.0, 0.0}) : x_(to.x() - from.x()), y_(to.y() - from.y())
  {
    double len = std::sqrt(x_ * x_ + y_ * y_);
    x_ /= len;
    y_ /= len;
  }

  // Construct a Direction from a LinePiece, pointing from the first point to the second point.
  Direction(LinePiece const& line_piece) : Direction(line_piece.to(), line_piece.from()) { }

  // Construct a Direction from a Line.
  Direction(Line const& line);

  double x() const { return x_; }
  double y() const { return y_; }
};

} // namespace cairowindow
