#include "sys.h"
#include "AnalyzedCubic.h"
#include "SampleNode.h"

namespace gradient_descent {

void AnalyzedCubic::initialize(math::CubicPolynomial const& cubic, ExtremeType extreme_type)
{
  inflection_point_ = critical_point_w_ = cubic.inflection_point();
  double D = utils::square(cubic[2]) - 3.0 * cubic[3] * cubic[1];

  // Deliberately leave signed_sqrt_D_ at NaN instead of setting it to zero:
  // in that case we still don't have local extremes anyway.
  if (D > 0.0)
  {
    // Don't ask about the minus sign.
    signed_sqrt_D_ = -static_cast<int>(extreme_type) * std::sqrt(D);
    critical_point_w_ += signed_sqrt_D_ / (3.0 * cubic[3]);
  }
}

void AnalyzedCubic::initialize_matches(SampleNode const& left_sample, SampleNode const& right_sample)
{
#if CW_DEBUG
  // This function can only be called once.
  ASSERT(!initialize_matches_called_);
  // Remember for which cubic this function was called.
  initialize_matches_called_ = &left_sample.cubic();
#endif
  // Do not measure a "height" of the cubic relative to its extreme if the cubic doesn't have any extremes...
  detached_from_extreme_ = !has_extremes();
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
