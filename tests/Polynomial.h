#pragma once

#include "utils/almost_equal.h"
#include "utils/square.h"
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

  int get_zeroes(std::array<double, 2>& zeroes_out) const
  {
    // This can be at most a parabola.
    ASSERT(1 <= coefficients_.size() && coefficients_.size() <= 3);
    if (coefficients_.size() < 3)
    {
      if (coefficients_.size() < 2)
        return 0;
      zeroes_out[0] = -coefficients_[0] / coefficients_[1];
      return 1;
    }

    double D = utils::square(coefficients_[1]) - 4.0 * coefficients_[2] * coefficients_[0];
    if (D < 0.0)
      return 0;
    double delta = std::sqrt(D) / std::abs(2.0 * coefficients_[2]);
    double avg = -coefficients_[1] / (2.0 * coefficients_[2]);
    zeroes_out[0] = avg - delta;
    zeroes_out[1] = avg + delta;
    return utils::almost_equal(zeroes_out[0], zeroes_out[1], 1e-6) ? 1 : 2;
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
