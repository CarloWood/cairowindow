#pragma once

#include "CS.h"
#include "math/TranslationVector.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include "utils/to_string.h"
#endif

namespace cairowindow::cs {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

template<math::CS cs> //FIXME: remove 'math::' once this becomes namespace math::cs
class Size
{
 private:
  double width_;
  double height_;

 public:
  constexpr Size(double width, double height) : width_(width), height_(height) { }

  double width() const { return width_; }
  double height() const { return height_; }

  // Implicit converstion to math::TranslationVector<cs>, so we can pass a Size<cs> to math::Transform<>::translate.
  operator math::TranslationVector<cs>() const
  {
    // This is OK because width_ and height are in the cs coordinate system.
    return math::TranslationVector<cs>::create_from_cs_values(width_, height_);
  }

  friend Size operator*(double scale, Size const& size)
  {
    return {scale * size.width(), scale * size.height()};
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":(" << width_ << ", " << height_ << ")";
  }
#endif
};

} // namespace cairowindow::cs
