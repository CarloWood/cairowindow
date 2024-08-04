#pragma once

#include "ExtremeType.h"
#include "../CubicPolynomial.h"
#include "utils/has_print_on.h"
#include <limits>

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

static_assert(std::numeric_limits<double>::has_quiet_NaN, "Sorry, you need a quiet NaN.");

class AnalyzedCubic
{
 private:
  double inflection_point_;
  double signed_sqrt_D_{std::numeric_limits<double>::quiet_NaN()};
  double critical_point_w_;             // If extreme, then the one passed to initialize.

 public:
  AnalyzedCubic() = default;

  void initialize(math::CubicPolynomial const& cubic, ExtremeType extreme_type_);

  bool has_extremes() const
  {
    // signed_sqrt_D_ is left as NaN when it is zero.
    return !std::isnan(signed_sqrt_D_);
  }

  double get_extreme() const
  {
    // Only call this function if has_extremes() returns true.
    ASSERT(has_extremes());
    return critical_point_w_;
  }

#ifdef CWDEBUG
  double get_other_extreme() const
  {
    return 2.0 * inflection_point_ - critical_point_w_;
  }
#endif

  double height(double w, double d) const
  {
    // Let g(w) = y = a + bw + cw^2 + dw^3.
    // Then dy/dw = b + 2cw + 3dw^2.
    // Let e be a root of the derivate: g'(e) = 0 (aka, g has an extreme in e).
    // Aka e = (-c +/- sqrt(c^2 - 3bd)) / (3d).
    // If we apply the coordinate transformation:
    //   w = w' + e, g'(w') = g(w) - g(e)
    // shifting the extreme to the orgin, then g' has an extreme in w' = 0 with value g'(0) = 0.
    // Working out g' gives:
    //   g'(w') = a + b(w'+e) + c(w'+e)^2 + d(w'+e)^3 - (a + be + ce^2 + de^3) =
    //          = a + bw' + be + cw'^2 + 2cw'e + ce^2 + dw'^3 + 3dw'^2 e + 3dw'e^2 + de^3 - a - be - ce^2 - de^3 =
    //          = (b + 2ce + 3de^2)w' + (c + 3de)w'^2 + dw'^3 =
    //          = (dw' +/- sqrt(c^2 - 3bd)) w'^2
    //
    // In the transformed coordinate system the height is given by
    //   |g'(w')| = |d w' +/- sqrt(c^2 - 3bd)| w'^2
    //
    double w_prime = w - critical_point_w_;
    double h = (d * w_prime + signed_sqrt_D_) * utils::square(w_prime);
    return std::abs(h);
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{inflection_point:" << inflection_point_ <<
         ", signed_sqrt_D:" << signed_sqrt_D_ <<
         ", critical_point_w:" << critical_point_w_ << '}';
  }
#endif
};

} // namespace gradient_descent
