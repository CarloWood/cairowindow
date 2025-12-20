#pragma once

#include "Direction.h"
#include "math/Line.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include "utils/to_string.h"
#endif

namespace cairowindow::cs {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<CS cs>
class Line;

template<CS cs>
struct LineTypes
{
  static constexpr int n = 2;
  using scalar_type = double;
  using point_type = Point<cs>;
  using direction_type = Direction<cs>;
  using derived_type = Line<cs>;
};

template<CS cs>
class Line : public math::LineOps<LineTypes<cs>>
{
 private:
  math::Line<2> raw_;

 public:
  // Construct an uninitialized line.
  Line() = default;

  // Construct a line through point with direction.
  Line(Point<cs> const& point, Direction<cs> const& direction) : raw_(point.raw(), direction.raw()) { }

  // Explicit conversion from math::Line.
  explicit Line(math::Line<2> const& raw) : raw_(raw) { }

  math::Line<2>& raw() { return raw_; }
  math::Line<2> const& raw() const { return raw_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":" << raw_;
  }
#endif
};

} // namespace cairowindow::cs
