#pragma once

#include "Point.h"
#include "math/Direction.h"

namespace cairowindow::cs {

template<CS cs>
class Direction;

template<CS cs>
class Vector;

template<CS cs>
struct DirectionTypes
{
  static constexpr int n = 2;
  using scalar_type = double;
  using vector_type = Vector<cs>;
  using derived_type = Direction<cs>;
};

template<CS cs>
class Direction : public math::DirectionOps<DirectionTypes<cs>>
{
 private:
  math::Direction<2> raw_;

 public:
  // Construct an undefined Direction.
  Direction() = default;

  // Construct a Direction that points in the direction theta (in radians): an angle with the positive x-axis.
  explicit Direction(double theta) : raw_(theta) { }

  // Construct a Direction from two points.
  Direction(Point<cs> const& from, Point<cs> const& to) : raw_(from.raw(), to.raw()) { }

  // If only one point is given, the direction is from the origin to that point.
  explicit Direction(Point<cs> const& to) : raw_(to.raw()) { }

  // Explicit conversion from math::Direction.
  explicit Direction(math::Direction<2> const& raw) : raw_(raw) { }

  math::Direction<2>& raw() { return raw_; }
  math::Direction<2> const& raw() const { return raw_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{d_:" << raw_ << '}';
  }
#endif
};

} // namespace cairowindow::cs
