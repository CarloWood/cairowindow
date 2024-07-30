#pragma once

#include "CriticalPointType.h"
#include "Sample.h"
#include "utils/square.h"
#include "debug.h"

namespace gradient_descent {

class Scale
{
 private:
  double value_{};                      // Cache of value(). An indication of what changes to w are significant,
                                        // relative to the critical point. This value is always larger than zero.
  CriticalPointType type_{CriticalPointType::none};     // Whether distances are measured relative to the minimum, maximum or
                                        // the inflection point (if the cubic has no extremes).
  double critical_point_w_;             // Cached value of the x-coordinate of the critical point.

  Sample const* left_edge_;             // Pointer to left-most Sample that is known to still match this cubic (within 10%).
  Sample const* right_edge_;            // Pointer to right-most Sample that is known to still match this cubic (within 10%).

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
  Scale(double value, CriticalPointType type, double critical_point_w, Sample const* left_edge, Sample const* right_edge) :
    value_(value), type_(type), critical_point_w_(critical_point_w), left_edge_(left_edge), right_edge_(right_edge) { }
  Scale& operator=(Scale const&) = default;

  // Accessors.
  double value() const { return value_; }
  CriticalPointType type() const { return type_; }
  double critical_point_w() const { return critical_point_w_; }
  double left_edge_w() const { return left_edge_->w(); }
  double right_edge_w() const { return right_edge_->w(); }

  void set(CriticalPointType type, double critical_point_w, Sample const* left_edge, Sample const* right_edge, double inflection_point_w)
  {
    type_ = type;
    critical_point_w_ = critical_point_w;
    left_edge_ = left_edge;
    right_edge_ = right_edge;
    value_ = calculate_value(inflection_point_w);
  }
};

} // namespace gradient_descent
