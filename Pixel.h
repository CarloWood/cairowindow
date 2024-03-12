#pragma once

#include <cmath>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace cairowindow {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

// Coordinates of a point in pixels (relative to the top-left corner of the Window).
class Pixel
{
 private:
  double x_;
  double y_;

 public:
  Pixel(int x, int y) : x_(x), y_(y) { }
  explicit Pixel(double x, double y) : x_(x), y_(y) { }

  double x() const { return x_; }
  double y() const { return y_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

} // namespace cairowindow
