#pragma once

#include "utils/almost_equal.h"
#include "utils/square.h"
#include "utils/macros.h"
#include <vector>
#include <array>
#include <ranges>
#include <cmath>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include <string>
#endif

// This should become part of machine-learning in the end.
namespace math {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Polynomial
{
 protected:
  std::vector<double> coefficients_;
#ifdef CWDEBUG
  std::string symbol_name_;
#endif

 public:
  // Create a polynomial of at most degree `number_of_coefficients - 1`, starting
  // with all coefficients set to zero.
  Polynomial(int number_of_coefficients COMMA_CWDEBUG_ONLY(std::string const& symbol_name)) :
    coefficients_(number_of_coefficients) COMMA_CWDEBUG_ONLY(symbol_name_(symbol_name)) { }

  double operator[](int i) const { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }
  double& operator[](int i) { ASSERT(0 <= i && i < coefficients_.size()); return coefficients_[i]; }

  double operator()(double w) const
  {
    double result = 0.0;
    for (double coefficient : std::ranges::reverse_view(coefficients_))
      result = w * result + coefficient;
    return result;
  };

  Polynomial derivative() const
  {
    Polynomial result(coefficients_.size() - 1 COMMA_CWDEBUG_ONLY(symbol_name_));
    for (int i = 1; i < coefficients_.size(); ++i)
      result[i - 1] = coefficients_[i] * i;
    return result;
  }

  // Return the division of this Polynomial by the factor (w - z).
  Polynomial long_division(double z, double& remainder) const
  {
    // f(w) = 3 * w^3 +      5 * w^2 - 4 * w + 10.
    //        3 * w^3 + (-2)*3 * w^2
    //      - ------------------------------------
    //                      11 * w^2 -     4 * w + 10.
    //                      11 * w^2 + (-2)*11 w
    //                    - --------------------------
    //                                    18 * w +     10.
    //                                    18 * w + (-2)18
    //                                  - ---------------
    //                                                 46
    // Divide by (w - 2)
    // 3 * w^2 + 11 * w + 18

    // NOTICE        : 10 + -4 w + 5 w^2 + 3 w^3
    // (w - 2)(3 w^2 + 11 w + 18) = 3 w^3 + 5 w^2 - 4 w - 36

    if (coefficients_.size() < 2)
    {
      ASSERT(coefficients_.size() == 1);
      remainder = coefficients_[0];
      return {1 COMMA_CWDEBUG_ONLY(symbol_name_)};
    }
    Polynomial result(coefficients_.size() - 1 COMMA_CWDEBUG_ONLY(symbol_name_));
    result[coefficients_.size() - 2] = coefficients_[coefficients_.size() - 1];
    for (int i  = coefficients_.size() - 2; i > 0; --i)
      result[i - 1] = coefficients_[i] + z * result[i];
    remainder = coefficients_[0] + z * result[0];
    return result;
  }

  int get_roots(std::array<double, 2>& roots_out) const
  {
    // This can be at most a parabola.
    ASSERT(1 <= coefficients_.size() && coefficients_.size() <= 3);
    if (coefficients_.size() < 3 || coefficients_[2] == 0.0)
    {
      if (coefficients_.size() < 2)
        return 0;
      roots_out[0] = -coefficients_[0] / coefficients_[1];
      return std::isfinite(roots_out[0]) ? 1 : 0;
    }

    double const D = utils::square(coefficients_[1]) - 4.0 * coefficients_[2] * coefficients_[0];
    if (D < 0.0)
      return 0;
    // Use a sqrt with the same sign as coefficients_[1];
    double const signed_sqrt_D = std::copysign(std::sqrt(D), coefficients_[1]);

    // Calculate the root closest to zero.
    roots_out[0] = -2.0 * coefficients_[0] / (coefficients_[1] + signed_sqrt_D);

    if (AI_UNLIKELY(std::isnan(roots_out[0])))
    {
      // This means we must have divided by zero, which means that both, coefficients_[1] as well as sqrtD, must be zero.
      // The latter means that coefficients_[0] is zero (coefficients_[2] was already checked not to be zero).
      // Therefore we have: f(x) = c x^2 with one root at x=0.
      roots_out[0] = 0.0;
      return 1;
    }

    // Calculate the root further away from zero.
    roots_out[1] = -0.5 * (coefficients_[1] + signed_sqrt_D) / coefficients_[2];

    // The second one is larger in absolute value.
    ASSERT(std::abs(roots_out[1]) > std::abs(roots_out[0]));

    return 2;
  }

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
          os << ' ' << symbol_name_;
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
