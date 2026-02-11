#pragma once

#include "CS.h"
#include "math/cs/Point.h"
#include "math/cs/Size.h"
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace cairowindow::cs {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<CS cs>
class Rectangle
{
 private:
  math::cs::Point<cs> top_left_;
  math::cs::Size<cs> size_;

 public:
  Rectangle() : top_left_{}, size_{0.0, 0.0} { }
  Rectangle(double offset_x, double offset_y, double width, double height) : top_left_{offset_x, offset_y}, size_{width, height}
  {
    ASSERT(width >= 0.0 && height >= 0.0);
  }
  Rectangle(math::cs::Point<cs> const& top_left, math::cs::Size<cs> const& size) : top_left_{top_left}, size_{size} { }

  math::cs::Point<cs> const& top_left() const { return top_left_; }
  math::cs::Size<cs> const& size() const { return size_; }

  double offset_x() const { return top_left_.x(); }
  double offset_y() const { return top_left_.y(); }
  double width() const { return size_.width(); }
  double height() const { return size_.height(); }

  bool is_defined() const { return size_.width() > 0.0 && size_.height() > 0.0; }

  double area() const { ASSERT(is_defined()); return size_.width() * size_.height(); }

  bool contains(int x, int y) const
  {
    return top_left_.x() <= x && x < (top_left_.x() + size_.width()) && top_left_.y() <= y && y < top_left_.y() + size_.height();
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":" << '{' << top_left_.x() << ", " << top_left_.y() << ", " << size_.width() << ", " << size_.height() << '}';
  }
#endif
};

} // namespace cairowindow::cs
