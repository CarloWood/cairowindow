#pragma once

#include <cmath>

#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include <iostream>
#endif

namespace cairowindow {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Range
{
 private:
  double min_;
  double max_;

 public:
  Range() : min_{0.0}, max_{1.0} { }
  Range(double min, double max) : min_(min), max_(max) { }

  // Increase the range so that both min_ and max_ are a multiple of spacing.
  void round_to(double spacing)
  {
    double const epsilon = 1e-4;
    min_ = spacing * std::floor(min_ / spacing + epsilon);
    max_ = spacing * std::ceil(max_ / spacing - epsilon);
  }

  // Accessors.
  double min() const { return min_; }
  double max() const { return max_; }
  double size() const { return max_ - min_; }
  double center() const { return 0.5 * (min_ + max_); }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << '[' << min_ << ", " << max_ << ']';
  }
#endif
};

} // namespace cairowindow
