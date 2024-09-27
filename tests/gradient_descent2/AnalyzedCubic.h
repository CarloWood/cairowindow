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

class SampleNode;

class AnalyzedCubic
{
 private:
  double inflection_point_;
  union {
    double signed_sqrt_D_{std::numeric_limits<double>::quiet_NaN()};
    double vertical_scale_;
  };
  double critical_point_w_;             // If extreme, then the one passed to initialize.
  bool detached_from_extreme_{false};   // Set if matches shouldn't use height. Only valid after calling initialize_matches.
#if CW_DEBUG
  bool initialize_matches_called_{false};
  math::CubicPolynomial const* debug_cubic_{nullptr};
#endif

 public:
  AnalyzedCubic() = default;

  void initialize(math::CubicPolynomial const& cubic, ExtremeType extreme_type_);
  void initialize_matches(SampleNode const& left_sample, SampleNode const& right_sample);

  bool has_extrema() const
  {
    // If this is set then signed_sqrt_D_ can't be used anymore (vertical_scale_ has been set).
    ASSERT(!detached_from_extreme_);
    // signed_sqrt_D_ is left as NaN when it is zero.
    return !std::isnan(signed_sqrt_D_);
  }

  double get_extreme() const
  {
    // Only call this function if has_extrema() returns true.
    ASSERT(has_extrema());
    return critical_point_w_;
  }

  double inflection_point() const
  {
    return inflection_point_;
  }

  double get_other_extreme() const
  {
    return 2.0 * inflection_point_ - critical_point_w_;
  }

  double height(double w, double d) const
  {
    // Only call this function if has_extrema() returned true.
    ASSERT(has_extrema());
    // Don't call height in this case.
    ASSERT(std::isnormal(signed_sqrt_D_));
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

  bool matches(double w, double Lw, math::CubicPolynomial const& g)
  {
    // If the vertical_scale_ is +inf then (w, Lw) always matches.
    if (!std::isfinite(vertical_scale_))
      return true;
    // Call initialize_matches before calling matches.
    ASSERT(initialize_matches_called_ && debug_cubic_ == &g);
    double const vertical_scale = detached_from_extreme_ ? vertical_scale_ : height(w, g[3]);
    // We want to return true if the point (w, Lw), deviates from the value according to the cubic g(w)
    // less than 10% of the (vertical) distance to the extreme: |g(w) - g(e)|.
    // In other words: |Lw - g(w)| <= 0.1 |g(w) - g(e)|
    return std::abs(Lw - g(w)) <= 0.1 * vertical_scale;
  }

#ifdef CWDEBUG
  // Return a pointer to the cubic that initalize was called with.
  math::CubicPolynomial const* debug_cubic() const { return debug_cubic_; }

  void print_on(std::ostream& os) const
  {
    os << "{inflection_point:" << inflection_point_ <<
         ", signed_sqrt_D:" << signed_sqrt_D_ <<
         ", critical_point_w:" << critical_point_w_ << '}';
  }
#endif
};

} // namespace gradient_descent
