#include "sys.h"
#include "AnalyzedCubic.h"

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

} // namespace gradient_descent
