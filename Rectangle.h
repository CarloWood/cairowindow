#pragma once

#include <cairo/cairo.h>
#include <memory>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace cairowindow {
#ifdef CWDEBUG
// This class defines a print_on method.
using utils::has_print_on::operator<<;
#endif

//FIXME: shouldn't this be coordinate system aware?
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
    os << '{' << geometry_.x << ", " << geometry_.y << ", " << geometry_.width << ", " << geometry_.height << '}';
  }
#endif
};

namespace draw {
class Rectangle;
} // namespace draw

namespace plot {
class Plot;

class Rectangle : public cairowindow::Rectangle
{
 public:
  explicit Rectangle(cairowindow::Rectangle const& rectangle) : cairowindow::Rectangle(rectangle) { }
  using cairowindow::Rectangle::Rectangle;

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Rectangle> draw_object_;
};

} // namespace plot
} // namespace cairowindow
