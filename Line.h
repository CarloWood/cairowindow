#pragma once

#include "math/Line.h"
#include <memory>

namespace cairowindow {
using Line = math::Line<2>;

namespace draw {
class Line;
} // namespace draw

namespace plot {
class Plot;

//--------------------------------------------------------------------------
// Line

class Line : public math::Line<2>
{
 public:
  using math::Line<2>::Line;
  Line(math::Line<2> const& line) : math::Line<2>(line) { }

  void reset()
  {
    draw_object_.reset();
  }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

} // namespace plot
} // namespace cairowindow
