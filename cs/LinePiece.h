#pragma once

#include "Direction.h"
#include "math/LinePiece.h"

namespace cairowindow::cs {

template<CS cs>
class LinePiece;

template<CS cs>
struct LinePieceTypes
{
  static constexpr int n = 2;
  using scalar_type = double;
  using point_type = Point<cs>;
  using direction_type = Direction<cs>;
  using derived_type = LinePiece<cs>;
};

template<CS cs>
class LinePiece : public math::LinePieceOps<LinePieceTypes<cs>>
{
 private:
  math::LinePiece<2> raw_;

 public:
  // Construct an uninitialized line.
  LinePiece() = default;

  // Construct a LinePiece from two points.
  LinePiece(Point<cs> const& from, Point<cs> const& to) : raw_(from.raw(), to.raw()) { }

  // Explicit conversion from math::LinePiece.
  explicit LinePiece(math::LinePiece<2> const& raw) : raw_(raw) { }

  math::LinePiece<2>& raw() { return raw_; }
  math::LinePiece<2> const& raw() const { return raw_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":" << raw_;
  }
#endif
};

} // namespace cairowindow::cs
