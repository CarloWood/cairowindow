#pragma once

#include <array>
#include <cmath>
#include "debug.h"
#include "utils/square.h"
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

  void initialize(double x0, double y0, double dxdy0, double x1, double y1, double dxdy1)
  {
    double delta_x = x0 - x1;
    double delta_x_squared = utils::square(delta_x);
    double delta_y = y0 - y1;
    double delta_dxdy = dxdy0 - dxdy1;
    double sum_x = x0 + x1;
    double sum_dxdy = dxdy0 + dxdy1;

    // The theory of this approach is described here:
    // https://math.stackexchange.com/a/4926903/489074
    double d = (-2.0 * delta_y + delta_x * sum_dxdy) / (delta_x_squared * delta_x);
    double c = delta_dxdy / (2.0 * delta_x) - 1.5 * sum_x * d;
    double b = (x0 * dxdy1 - x1 * dxdy0) / delta_x + 3.0 * x0 * x1 * d;

    coefficients_[3] = d;
    coefficients_[2] = c;
    coefficients_[1] = b;
    coefficients_[0] = 0.0;
    coefficients_[0] = y0 - operator()(x0);
  }

  int get_extremes(std::array<double, 2>& extremes_out)
  {
    DoutEntering(dc::notice, "CubicPolynomial::get_extremes()");

    double D = utils::square(coefficients_[2]) - 3.0 * coefficients_[1] * coefficients_[3];
    Dout(dc::notice, "D = " << D);

    if (D < 0.0)
      return 0;

    int index_minimum = (coefficients_[3] > 0.0) ? 1 : 0;

    double sqrt_D = std::sqrt(D);
    extremes_out[index_minimum] = (-coefficients_[2] + sqrt_D) / (3.0 * coefficients_[3]);
    extremes_out[1 - index_minimum] = (-coefficients_[2] - sqrt_D) / (3.0 * coefficients_[3]);

    Dout(dc::notice, "extremes_out = " << extremes_out);

    return (D == 0.0) ? 1 : 2;
  }

  // Evaluation.
  double operator()(double x) const
  {
    return coefficients_[0] + (coefficients_[1] + (coefficients_[2] + coefficients_[3] * x) * x) * x;
  };

  // Evaluate derivative.
  double derivative(double x) const
  {
    return coefficients_[1] + (2.0 * coefficients_[2] + 3.0 * coefficients_[3] * x) * x;
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
