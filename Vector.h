#pragma once

#include "math/Vector.h"
#include "Pixel.h"

namespace cairowindow {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Vector : public math::Vector
{
 public:
  using math::Vector::Vector;
  Vector(math::Vector const& v) : math::Vector(v) { }

  // Convert the vector to a Pixel.
  Pixel as_pixel() const { return Pixel{x_, y_}; }
};

} // namespace cairowindow
