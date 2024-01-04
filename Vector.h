#pragma once

#include "LinePiece.h"
#include <cmath>

namespace cairowindow {

class Vector
{
 private:
  double x_;
  double y_;

 public:
  // Construct a vector from its x,y coordinates.
  Vector(double x, double y) : x_(x), y_(y) { }

  // Construct a Vector that points in direction and has length.
  // Also used for automatic conversion from a Direction to a Vector.
  Vector(Direction direction, double length = 1.0) : x_(direction.x() * length), y_(direction.y() * length) { }

  // Construct a Vector from two points. If the second point is not given it defaults to the origin.
  // The direction is from the second argument (or origin) to the first argument.
  Vector(Point const& from, Point const& to) : x_(to.x() - from.x()), y_(to.y() - from.y()) { }
  explicit Vector(Point const& to) : x_(to.x()), y_(to.y()) { }

  // Construct a Vector from a LinePiece, pointing from the first point to the second point.
  Vector(LinePiece const& line_piece) : Vector(line_piece.from(), line_piece.to()) { }

  double x() const { return x_; }
  double y() const { return y_; }

  // Return dot product with d2.
  double dot(Vector const& v2) const { return x_ * v2.x_ + y_ * v2.y_; }

  // Construct a Direction from this vector.
  Direction direction() const { return Point{x_, y_}; }

  // Return the length of the vector.
  double length() const { return std::sqrt(x_ * x_ + y_ * y_); }

  // Convert the vector to a Point.
  Point point() const { return {x_, y_}; }

 public:
  // Return the vector rotated 90 degrees counter-clockwise.
  Vector rotate_90_degrees() const { return { -y_, x_ }; }

  // Return the vector rotated 180 degrees.
  Vector rotate_180_degrees() const { return { -x_, -y_ }; }

  // Return the vector rotated 270 degrees.
  Vector rotate_270_degrees() const { return { y_, -x_ }; }
};

inline Vector operator*(double length, Vector const& v2)
{
  return {length * v2.x(), length * v2.y()};
}

inline Point operator+(Point const& point, Vector const& v2)
{
  return {point.x() + v2.x(), point.y() + v2.y()};
}

inline Vector operator+(Vector const& v1, Vector const& v2)
{
  return {v1.x() + v2.x(), v1.y() + v2.y()};
}

} // namespace cairowindow
