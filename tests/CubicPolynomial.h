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
 protected:
  std::array<double, 4> coefficients_{};

 public:
  // Create a zero polynomial.
  CubicPolynomial() { }
  // Create a polynomial  a + b x + c x^2 + d x^3.
  CubicPolynomial(double a, double b, double c, double d) : coefficients_{{a, b, c, d}} { }

  void initialize(double x0, double y0, double dxdy0, double x1, double y1, double dxdy1)
  {
    double delta_x_inverse = 1.0 / (x0 - x1);

    // The theory of this approach is described here:
    // https://math.stackexchange.com/a/4926903/489074
    double d = (dxdy0 + dxdy1 - 2.0 * (y0 - y1) * delta_x_inverse) * (delta_x_inverse * delta_x_inverse);
    double c = 0.5 * ((dxdy0 - dxdy1) * delta_x_inverse - 3.0 * d * (x0 + x1));
    double b = (x0 * dxdy1 - x1 * dxdy0) * delta_x_inverse + 3.0 * d * x0 * x1;

    coefficients_[3] = d;
    coefficients_[2] = c;
    coefficients_[1] = b;
    coefficients_[0] = 0.0;
    coefficients_[0] = y0 - operator()(x0);
  }

  int get_extrema(std::array<double, 2>& extrema_out, bool left_most_first = true) const
  {
    DoutEntering(dc::notice, "CubicPolynomial::get_extrema()");

    double D = utils::square(coefficients_[2]) - 3.0 * coefficients_[1] * coefficients_[3];
    Dout(dc::notice, "D = " << D);

    if (D < 0.0)
    {
      // Write the inflection point to index 0.
      extrema_out[0] = -coefficients_[2] / (3.0 * coefficients_[3]);
      return 0;
    }

    // Put the left-most extreme in index 0 iff left_most_first is true.
    // Otherwise, put the minimum in index 0.
    int index_minimum = (left_most_first && coefficients_[3] > 0.0) ? 1 : 0;

    double sqrt_D = std::sqrt(D);
    extrema_out[index_minimum] = (-coefficients_[2] + sqrt_D) / (3.0 * coefficients_[3]);
    extrema_out[1 - index_minimum] = (-coefficients_[2] - sqrt_D) / (3.0 * coefficients_[3]);

    Dout(dc::notice, "extrema_out = " << std::setprecision(std::numeric_limits<double>::digits10) << extrema_out);

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

  // Evaluate second derivative.
  double second_derivative(double x) const
  {
    return 2.0 * coefficients_[2] + 6.0 * coefficients_[3] * x;
  }

  double inflection_point() const
  {
    return -coefficients_[2] / (3.0 * coefficients_[3]);
  }

  // Access coefficients.
  double operator[](int i) const { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }
  double& operator[](int i) { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }

#if CW_DEBUG
  // Return true if one was assigned from the other.
  friend bool operator==(CubicPolynomial const& lhs, CubicPolynomial const& rhs)
  {
    return lhs.coefficients_ == rhs.coefficients_;
  }
#endif

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
