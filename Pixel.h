#pragma once

#include "utils/has_print_on.h"
#include <cmath>

namespace cairowindow {
using utils::has_print_on::operator<<;

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
