#pragma once

#include "LinePiece.h"
#include <cmath>
#include "debug.h"

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

  // Return dot product with d2.
  double dot(Direction const& d2) const { return x_ * d2.x_ + y_ * d2.y_; }

  double as_angle() const { return std::atan2(y_, x_); }

 protected:
  // For normal() and inverse().
  constexpr Direction(double x, double y) : x_(x), y_(y) { };

 public:
  // Return the direction rotated 90 degrees counter-clockwise.
  Direction normal() const { return { -y_, x_ }; }

  // Return the direction rotated 180 degrees.
  Direction inverse() const { return { -x_, -y_ }; }

  static Direction const up;
  static Direction const down;
  static Direction const left;
  static Direction const right;
};

} // namespace cairowindow
