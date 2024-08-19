#pragma once

#include "CriticalPointType.h"
#include "Sample.h"
#include "AnalyzedCubic.h"
#include "HorizontalDirection.h"
#include "../CubicPolynomial.h"
#include "utils/square.h"
#include "debug.h"

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Scale
{
 private:
  double value_{};                      // Cache of value(). An indication of what changes to w are significant,
                                        // relative to the critical point. This value is always larger than zero.
  CriticalPointType type_{CriticalPointType::none};     // Whether distances are measured relative to the minimum, maximum or
                                        // the inflection point (if the cubic has no extremes).
  double critical_point_w_;             // Cached value of the x-coordinate of the critical point.

  double left_edge_w_;                  // The w value of left-most Sample that is known to still match this cubic (within 10%).
  double right_edge_w_;                 // The w value of right-most Sample that is known to still match this cubic (within 10%).

 public:
  static constexpr int Li = 0;
  static constexpr int Ci = 1;
  static constexpr int Ri = 2;
  static constexpr int Ii = 3;
  using LCRI_type = std::array<double, 4>;      // An array with the L, C, R and I values.

  static constexpr unsigned char const classify_map[13] = {4, 4, 1, 4, 3, 2, 0, 3, 4, 1, 4, 4};
  static int classify(LCRI_type const& LCRI)
  {
    return classify_map[(LCRI[Ii] < LCRI[Li]) + (LCRI[Ii] < LCRI[Ri]) +
      4 * ((LCRI[Ci] < LCRI[Li]) + (LCRI[Ci] < LCRI[Ri])) + (LCRI[Ii] < LCRI[Ci])];
  }

 private:
  // Return the weighted average distance of L and R to C.
  static double weighted_average(LCRI_type const& LCRI)
  {
    ASSERT(LCRI[Li] < LCRI[Ri]);        // In fact, LCRI[Li] < LCRI[Ci] < LCRI[Ri].
    DoutEntering(dc::notice, "Scale::weighted_average(" << LCRI[Li] << ", " << LCRI[Ci] << ", " << LCRI[Ri] << ")");
    return (utils::square(LCRI[Ci] - LCRI[Li]) + utils::square(LCRI[Ri] - LCRI[Ci])) / (LCRI[Ri] - LCRI[Li]);
  }

  double calculate_value(double const inflection_point_w) const;

 public:
  Scale() = default;
  Scale(Scale const&) = default;
  Scale(double value, CriticalPointType type, double critical_point_w, double left_edge_w, double right_edge_w) :
    value_(value), type_(type), critical_point_w_(critical_point_w), left_edge_w_(left_edge_w), right_edge_w_(right_edge_w) { }
  Scale& operator=(Scale const&) = default;

  // Accessors.
  double value() const { return value_; }
  CriticalPointType type() const { return type_; }
  double critical_point_w() const { return critical_point_w_; }
  double left_edge_w() const { return left_edge_w_; }
  double right_edge_w() const { return right_edge_w_; }

  void set(CriticalPointType type, double critical_point_w, double left_edge_w, double right_edge_w, double inflection_point_w)
  {
    ASSERT((type == CriticalPointType::inflection_point) == (critical_point_w == inflection_point_w));
    type_ = type;
    critical_point_w_ = critical_point_w;
    left_edge_w_ = left_edge_w;
    right_edge_w_ = right_edge_w;
    value_ = calculate_value(inflection_point_w);
  }

  // Return a directional scale: minus the scale if direction is `left` and plus the scale if direction is `right`.
  double step(HorizontalDirection& hdirection) const
  {
    DoutEntering(dc::notice, "Scale::step(" << hdirection << ")");
    ASSERT(type_ != CriticalPointType::none);
    // If hdirection is undecided at this point then make a step in the direction from the sample (part of scale)
    // that is the furthest away from the critical point, towards the critical point.
    if (hdirection == HorizontalDirection::undecided)
    {
      double far_edge_w =
        std::abs(left_edge_w_ - critical_point_w_) > std::abs(right_edge_w_ - critical_point_w_) ? left_edge_w_ : right_edge_w_;
      hdirection = critical_point_w_ < far_edge_w ? HorizontalDirection::left : HorizontalDirection::right;
    }
    return static_cast<int>(hdirection) * value_;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << '{';
    if (type_ != CriticalPointType::none)
    {
      os << "value:" << value_ <<
          ", type:" << type_ <<
          ", cp:" << critical_point_w_ <<
          ", left_edge_w:" << left_edge_w_ <<
          ", right_edge_w:" << right_edge_w_;
    }
    else
      os << "<no scale>";
    os << '}';
  }
#endif
};

} // namespace gradient_descent
