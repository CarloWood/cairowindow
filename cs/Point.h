#pragma once

#include "Size.h"
#include "math/TranslationVector.h"
#include "math/Point.h"

namespace cairowindow::cs {

template<CS cs>
class Point;

template<CS cs>
class Vector;

template<CS cs>
class Direction;

template<CS cs>
struct PointTypes
{
  static constexpr int n = 2;
  using scalar_type = double;
  using vector_type = Vector<cs>;
  using direction_type = Direction<cs>;
  using derived_type = Point<cs>;
};

template<CS cs>
class Point : public math::PointOps<PointTypes<cs>>
{
 private:
  math::Point<2> raw_;

 public:
  // Construct the origin.
  Point() : raw_(0.0, 0.0) { }

  // Construct a Point from two coordinates.
  constexpr Point(double x, double y) : raw_(x, y) { }

  //---------------------------------------------------------------------------
  // Constructors from math::Point.
  Point(typename math::Point<2>::eigen_type const& point) : raw_(point) { }
  explicit Point(math::Point<2> const& point) : raw_(point) { }
  //---------------------------------------------------------------------------

  // Add a Size<cs>.
  Point operator+(Size<cs> const& size) const
  {
    return {raw_.x() + size.width(), raw_.y() + size.height()};
  }
  // These did get hidden because of the above operator.
  using math::PointOps<PointTypes<cs>>::operator+;

  // Accessors to the underlying math::Point (required for math::PointOps).
  math::Point<2>& raw() { return raw_; }
  math::Point<2> const& raw() const { return raw_; }

  // Implicit converstion to math::TranslationVector<cs>, so we can pass a Point<cs> to math::Transform<>::translate.
  operator math::TranslationVector<cs>() const
  {
    return math::TranslationVector<cs>::create_from_cs_values(raw_);
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":(" << raw_ << ")";
  }
#endif
};

} // namespace cairowindow::cs
