#pragma once

#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include <iostream>
#endif

namespace cairowindow {
#ifdef CWDEBUG
// This class defines a print_on method.
using utils::has_print_on::operator<<;
#endif

class Color
{
 private:
  double red_;
  double green_;
  double blue_;
  double alpha_;

 public:
  constexpr Color(double red, double green, double blue, double alpha = 1.0) : red_(red), green_(green), blue_(blue), alpha_(alpha) { }

  double red() const { return red_; }
  double green() const { return green_; }
  double blue() const { return blue_; }
  double alpha() const { return alpha_; }

  bool is_transparent() const { return alpha_ == 0.0; }
  bool is_opaque() const { return alpha_ == 1.0; }
  void set_opaque() { alpha_ = 1.0; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{red_:" << red_ << ", green_:" << green_ << ", blue_:" << blue_ << ", alpha_:" << alpha_ << '}';
  }
#endif
};

namespace color {

constexpr Color red{1.0, 0.0, 0.0};
constexpr Color green{0.0, 1.0, 0.0};
constexpr Color blue{0.0, 0.0, 1.0};
constexpr Color yellow{1.0, 1.0, 0.0};
constexpr Color magenta{1.0, 0.0, 1.0};
constexpr Color cyan{0.0, 1.0, 1.0};
constexpr Color black{0.0, 0.0, 0.0};
constexpr Color white{1.0, 1.0, 1.0};
constexpr Color gray{0.5, 0.5, 0.5};
constexpr Color orange{1.0, 0.647, 0.0};
constexpr Color brown{0.647, 0.165, 0.165};
constexpr Color pink{1.0, 0.753, 0.796};
constexpr Color violet{0.933, 0.51, 0.933};
constexpr Color indigo{0.294, 0.0, 0.51};
constexpr Color maroon{0.502, 0.0, 0.0};
constexpr Color turquoise{0.251, 0.878, 0.816};
constexpr Color tan{0.824, 0.706, 0.549};
constexpr Color gold{1.0, 0.843, 0.0};
constexpr Color beige{0.961, 0.961, 0.863};
constexpr Color lavender{0.902, 0.902, 0.980};
constexpr Color peach{1.0, 0.855, 0.725};
constexpr Color lime{0.0, 1.0, 0.0};
constexpr Color olive{0.502, 0.502, 0.0};
constexpr Color teal{0.0, 0.502, 0.502};
constexpr Color navy{0.0, 0.0, 0.502};
constexpr Color coral{1.0, 0.498, 0.314};
constexpr Color salmon{0.980, 0.502, 0.447};

constexpr Color transparent{0.0, 0.0, 0.0, 0.0};

} // namespace color

} // namespace cairowindow
