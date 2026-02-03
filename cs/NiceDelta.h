#pragma once

#include "Range.h"
#include "math/subsuper_string.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <string>

namespace cairowindow::cs {

template<CS cs>
class NiceDelta
{
 public:
  static constexpr int number_of_mantissa_values = 3;
  static constexpr std::array<int, number_of_mantissa_values> mantissa_values = { 1, 2, 5 };
  static constexpr int invalid_magic = -2;

 private:
  int mantissa_;
  int exponent_;
  int m_;

 private:
  // Construct NiceDelta from a given mantissa and exponent; this constructor accepts a mantissa in the range [-1, number_of_mantissa_values - 1].
  // Note that m_ is not initialized. This object may only be used to call value() on it.
  NiceDelta(int mantissa, int exponent) : mantissa_(mantissa), exponent_(exponent)
  {
    // Cannonicalize; turn a mantissa value of -1 into number_of_mantissa_values - 1 and decrement exponent_.
    if (mantissa_ == -1)
    {
      mantissa_ += number_of_mantissa_values;
      --exponent_;
    }
  }

  void increment()
  {
    mantissa_ = (mantissa_ + 1) % number_of_mantissa_values;
    if (mantissa_ == 0)
      ++exponent_;
  }

  static int calculate_m(Range<cs> const& range, double delta)
  {
    return static_cast<int>(std::floor(range.max() / delta) - std::ceil(range.min() / delta)) + 1;
  }

 public:
  // Construct an "invalid" NiceDelta.
  NiceDelta() : mantissa_(invalid_magic) { }

  NiceDelta(Range<cs> const& range)
  {
    // Dividing by 9 guarantees that m() will return a value less than 10,
    // which then is guaranteed not larger than the sought for value.
    double ideal_delta = range.size() / 9;

    // Initialize starting value for mantissa_ and exponent_.
    double ideal_exponent = std::log10(ideal_delta);
    exponent_ = static_cast<int>(std::floor(ideal_exponent));
    mantissa_ = 0;

    // Increment the current value until it is no longer less than ideal_delta.
    double current_value;
    for (;;)
    {
      current_value = value();

      // Exit once we found the first value that is greater or equal ideal_delta.
      if (current_value >= ideal_delta)
        break;

      // Go to the next allowed value.
      increment();
    }

    // Calculate m for the current value.
    m_ = calculate_m(range, current_value);

    // A value of 6 or larger is the correct value: a smaller delta would lead to an m of more than 10.
    if (m_ <= 5)
    {
      // Try the next smaller value of the current NiceDelta.
      NiceDelta next_smaller_delta(mantissa_ - 1, exponent_);
      int next_m = calculate_m(range, next_smaller_delta.value());

      if (m_ < 5 || next_m <= 10)
      {
        exponent_ = next_smaller_delta.exponent_;
        mantissa_ = next_smaller_delta.mantissa_;
        m_ = next_m;
      }
    }
  }

  bool is_invalid() const
  {
    return mantissa_ == invalid_magic;
  }

  double value() const
  {
    return mantissa_values[mantissa_] * std::pow(10.0, exponent_);
  }

  int m() const
  {
    return m_;
  }

  std::string label(int k) const
  {
    ASSERT(!is_invalid());

    double const delta_value = value();
    double const tick_value = k * delta_value;

    if (tick_value == 0.0)
      return "0";

    if (exponent_ <= -4 || exponent_ > 4)
    {
      long long const scaled_integer = static_cast<long long>(k) * mantissa_values[mantissa_];
      if (scaled_integer == 0)
        return "0";
      std::ostringstream oss;
      oss << scaled_integer << 'e' << exponent_;
      return oss.str();
    }

    std::ostringstream oss;
    oss.setf(std::ios::fixed);
    oss << std::setprecision(std::max(0, -exponent_)) << tick_value;
    std::string result = oss.str();
    if (!result.empty() && result.front() == '-')
    {
      auto const non_zero = std::find_if(result.begin() + 1, result.end(), [](char c) { return c != '0' && c != '.'; });
      if (non_zero == result.end())
        result.erase(result.begin());
    }
    return result;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    if (is_invalid())
      os << "<invalid>";
    else
      os << utils::to_string(cs) << ":{" << mantissa_values[mantissa_] << "·10" << math::to_superscript(exponent_) << "; m:" << m_ << "}";
  }
#endif
};

} // namespace cairowindow::cs
