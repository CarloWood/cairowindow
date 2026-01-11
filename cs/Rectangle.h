#pragma once

#include "Point.h"              // Point, Size
#include <cairo/cairo.h>
#include "debug.h"

namespace cairowindow::cs {

template<CS cs>
class Rectangle
{
 private:
  cairo_rectangle_t geometry_;

 public:
  Rectangle() : geometry_{} { }
  Rectangle(double offset_x, double offset_y, double width, double height) : geometry_{offset_x, offset_y, width, height}
  {
    ASSERT(width >= 0.0 && height >= 0.0);
  }

  double offset_x() const { return geometry_.x; }
  double offset_y() const { return geometry_.y; }
  double width() const { return geometry_.width; }
  double height() const { return geometry_.height; }

  bool is_defined() const { return geometry_.width > 0.0 && geometry_.height > 0.0; }

  double area() const { ASSERT(is_defined()); return geometry_.width * geometry_.height; }

  bool contains(int x, int y) const
  {
    return geometry_.x <= x && x < (geometry_.x + geometry_.width) && geometry_.y <= y && y < geometry_.y + geometry_.height;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":" << '{' << geometry_.x << ", " << geometry_.y << ", " << geometry_.width << ", " << geometry_.height << '}';
  }
#endif
};

} // namespace cairowindow::cs
