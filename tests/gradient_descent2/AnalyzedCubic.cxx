#include "sys.h"
#include "AnalyzedCubic.h"
#include "SampleNode.h"

namespace gradient_descent {

void AnalyzedCubic::initialize(math::CubicPolynomial const& cubic, ExtremeType extreme_type)
{
#if CW_DEBUG
  // Remember for which cubic this function was called.
  debug_cubic_ = &cubic;
#endif

  double half_sqrt_D;
  if (AI_UNLIKELY(cubic[3] == 0.0))
  {
    // We don't really have an inflection point.
    inflection_point_ = std::numeric_limits<double>::infinity();

    // In this case the cubic is a parabola:
    //
    //   A(w) = b w^2 + c w + d
    //
    // with derivative:
    //
    //   A'(w) = 2b w + c
    //
    // The second derivative is a constant:
    //
    //   A''(w) = 2b
    //
    // Hence whether or not r is a minimum or maximum depends on the sign of b.
    // If b < 0 then it is a maximum.

    if ((extreme_type == ExtremeType::maximum) == (cubic[2] < 0.0))
    {
      // The derivative has one root at -c / 2b.
      critical_point_w_ = -0.5 * cubic[1] / cubic[2];
      // See below, assuming cubic[3] = 0.
      half_sqrt_D = std::abs(cubic[2]);
    }
    else
    {
      // Pretend the critical_point_w_ is at plus infinity.
      critical_point_w_ = std::numeric_limits<double>::infinity();
      // We should never be using this to calculate a height...
      half_sqrt_D = std::numeric_limits<double>::infinity();
    }
  }
  else
  {
    inflection_point_ = cubic.inflection_point();

    double one_fourth_D = utils::square(cubic[2]) - 3.0 * cubic[1] * cubic[3];
    if (one_fourth_D <= 0.0)
    {
      // If the determinant is zero, then the cubic has no local extremes.
      return;
    }

    // Use a sqrt with the same sign as cubic[2];
    half_sqrt_D = std::sqrt(one_fourth_D);
    double const half_Q = std::copysign(half_sqrt_D, cubic[2]);

    // The roots of the derivative are:
    // x_0 = -c / (b + 0.5 * Q);
    // x_1 = (-b - 0.5 * Q) / (3 * a);
    // where
    // f''(x_0) = Q
    // f''(x_1) = -Q
    // Therefore if Q is positive then x_0 is the minimum and x_1 is the maximum
    // and if Q is negative then x_0 is the maximum and x_1 is the minimum.
    //
    // Note: if cubic[2] is zero (or close to zero due to floating point round of errors)
    // then the sign of half_Q is not well defined, but its absolute value is usually still
    // significant. However, in that case, x_0 == -x_1 and changing the sign of half_Q has
    // no real influence because that is exactly where we swap formula as well.
    if ((extreme_type == ExtremeType::maximum) == (half_Q < 0.0))
    {
      // Calculate the root closest to zero.
      critical_point_w_ = -cubic[1] / (cubic[2] + half_Q);
    }
    else
    {
      // Calculate the root further away from zero.
      critical_point_w_ = -(cubic[2] + half_Q) / (3.0 * cubic[3]);
    }
  }
  // Don't ask about the minus sign.
  signed_sqrt_D_ = -2 * static_cast<int>(extreme_type) * half_sqrt_D;
}

void AnalyzedCubic::initialize_matches(SampleNode const& left_sample, SampleNode const& right_sample)
{
#if CW_DEBUG
  // Call initialize before calling initialize_matches.
  ASSERT(debug_cubic_);
  // This function can only be called once.
  ASSERT(!initialize_matches_called_);
  initialize_matches_called_ = true;
  // This should match.
  ASSERT(debug_cubic_ == &left_sample.cubic());
#endif
  // Do not measure a "height" of the cubic relative to its extreme if the cubic doesn't have any extrema...
  detached_from_extreme_ = !has_extrema();
  if (detached_from_extreme_)
    vertical_scale_ = std::abs(left_sample.Lw() - right_sample.Lw());
  else
  {
    double const left_w = left_sample.w();
    double const right_w = right_sample.w();
    double const extreme_w = get_extreme();
    // but also not if the extreme does not fall within the two samples used for the cubic...
    if (extreme_w < left_w || right_w < extreme_w)
    {
      math::CubicPolynomial const& cubic = left_sample.cubic();
      double const left_Lw = left_sample.Lw();
      double const right_Lw = right_sample.Lw();
      double const extreme_Lw = cubic(extreme_w);
      double const center_Lw = 0.5 * (left_Lw + right_Lw);
      double const delta_Lw = std::abs(left_Lw - right_Lw);
      double const extreme_to_center_Lw = std::abs(extreme_Lw - center_Lw);
      // and the extreme is verticaly further away from the nearest sample then the vertical distance between the two samples.
      detached_from_extreme_ = extreme_to_center_Lw > 1.5 * delta_Lw;
      if (detached_from_extreme_)
        vertical_scale_ = delta_Lw;
    }
  }
}

} // namespace gradient_descent
