#pragma once

#include "CS.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include "utils/to_string.h"
#endif

namespace cairowindow::cs {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<CS cs>
class Size
{
 private:
  double width_;
  double height_;

 public:
  constexpr Size(double width, double height) : width_(width), height_(height) { }

  double width() const { return width_; }
  double height() const { return height_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":(" << width_ << ", " << height_ << ")";
  }
#endif
};

} // namespace cairowindow::cs
