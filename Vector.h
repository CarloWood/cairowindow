#pragma once

#include "math/Vector.h"
#include "Pixel.h"

namespace cairowindow {

class Vector : public math::Vector<2>
{
 public:
  using math::Vector<2>::Vector;
  Vector(math::Vector<2> const& v) : math::Vector<2>(v) { }

  // Convert the vector to a Pixel.
  Pixel as_pixel() const { return Pixel{x(), y()}; }
};

} // namespace cairowindow
