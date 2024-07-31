#include "sys.h"
#include "AnalyzedCubic.h"

namespace gradient_descent {

AnalyzedCubic::AnalyzedCubic(math::CubicPolynomial const& cubic) :
  math::CubicPolynomial(cubic), inflection_point_(inflection_point())
{
  double D = utils::square(coefficients_[2]) - 3.0 * coefficients_[3] * coefficients_[1];

  // Deliberately leave sqrt_D_div_3_ at NaN instead of setting it to zero:
  // in that case we still don't have local extremes anyway.
  if (D > 0.0)
    sqrt_D_div_3_ = std::sqrt(D) / (3.0 * coefficients_[3]);
}

} // namespace gradient_descent
