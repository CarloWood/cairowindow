#pragma once

#include "LinePiece.h"
#include "cairowindow/Pixel.h"
#include "math/Vector.h"

namespace cairowindow::cs {

template<CS cs>
class Vector;

template<CS cs>
struct VectorTypes
{
  static constexpr int n = 2;
  using scalar_type = double;
  using point_type = Point<cs>;
  using direction_type = Direction<cs>;
  using derived_type = Vector<cs>;
};

template<CS cs>
class Vector : public math::VectorOps<VectorTypes<cs>>
{
 private:
  math::Vector<2> raw_;

 public:
  // Construct an uninitialized Vector.
  Vector() = default;

  // Construct a standard basis vector e_i, or zero if i == -1.
  Vector(int i) : raw_(i) { }

  // Construct a vector from its two coordinates. The caller is responsible to make sure that the passed coordinates are from the right coordinate system.
  Vector(double x, double y) : raw_(x, y) { }

  // Construct a Vector that points in direction and has length.
  // Also used for automatic conversion from a Direction to a Vector.
  Vector(Direction<cs> direction, double length = 1.0) : raw_(direction.raw(), length) { }

  // Construct a Vector from two points. If the second point is not given it defaults to the origin.
  // The direction is from the second argument (or origin) to the first argument.
  Vector(Point<cs> const& from, Point<cs> const& to) : raw_(from.raw(), to.raw()) { }
  explicit Vector(Point<cs> const& to) : raw_(to.raw()) { }

  // Construct a Vector from a Line, pointing from the first point to the second point.
  explicit Vector(LinePiece<cs> const& line_piece) : raw_(line_piece.raw()) { }

  // Explicit conversion from math::Vector.
  explicit Vector(math::Vector<2> const& raw) : raw_(raw) { }

  // Convert the vector to a Pixel.
  Point<CS::pixels> as_pixel() const { return {this->x(), this->y()}; }

  math::Vector<2>& raw() { return raw_; }
  math::Vector<2> const& raw() const { return raw_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":" << raw_;
  }
#endif
};

} // namespace cairowindow::cs
