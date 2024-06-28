#pragma once

#include <array>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class CubicPolynomial
{
 private:
  std::array<double, 4> coefficients_{};

 public:
  // Create a zero polynomial.
  CubicPolynomial() { }
  // Create a polynomial  a + b x + c x^2 + d x^3.
  CubicPolynomial(double a, double b, double c, double d) : coefficients_{{a, b, c, d}} { }

  // Evaluation.
  double operator()(double w) const
  {
    return coefficients_[0] + (coefficients_[1] + (coefficients_[2] + coefficients_[3] * w) * w) * w;
  };

  // Evaluate derivative.
  double derivative(double w) const
  {
    return coefficients_[1] + (2.0 * coefficients_[2] + 3.0 * coefficients_[3] * w) * w;
  }

  // Access coefficients.
  double operator[](int i) const { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }
  double& operator[](int i) { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    bool first = true;
    int exponent = 0;
    for (double coefficient : coefficients_)
    {
      if (coefficient != 0.0)
      {
        if (first)
          os << coefficient;
        else if (coefficient > 0.0)
          os << " + " << coefficient;
        else
          os << " - " << -coefficient;
        if (exponent > 0)
        {
          os << " x";
          if (exponent > 1)
            os << '^' << exponent;
        }
        first = false;
      }
      ++exponent;
    }
  }
#endif
};

} // namespace math
