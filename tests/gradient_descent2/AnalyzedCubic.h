#pragma once

#include "ExtremeType.h"
#include "../CubicPolynomial.h"
#include <limits>

namespace gradient_descent {

static_assert(std::numeric_limits<double>::has_quiet_NaN, "Sorry, you need a quiet NaN.");

class AnalyzedCubic : public math::CubicPolynomial
{
 private:
  double inflection_point_;
  double sqrt_D_div_3_{std::numeric_limits<double>::quiet_NaN()};

 public:
  AnalyzedCubic(math::CubicPolynomial const& cubic);

  bool has_extremes() const
  {
    // sqrt_D_div_3_ is left as NaN when it is zero.
    return !std::isnan(sqrt_D_div_3_);
  }

  double get_extreme(ExtremeType extreme_type_) const
  {
    // Only call this function if has_extremes() returns true.
    ASSERT(has_extremes());

    return inflection_point_ - static_cast<int>(extreme_type_) * sqrt_D_div_3_;
  }
};

} // namespace gradient_descent
