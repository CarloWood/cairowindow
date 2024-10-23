#include "sys.h"
#include "AnalyzedCubic.h"
#include "SampleNode.h"

namespace gradient_descent {

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
  // signed_sqrt_D_ is used for something else if detached_from_extreme_ is set.
  double& vertical_scale_ = signed_sqrt_D_;
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
