#pragma once

#include "math/Line.h"
#include <memory>

namespace cairowindow {
using Line = math::Line;

namespace draw {
class Line;
} // namespace draw

namespace plot {
class Plot;

//--------------------------------------------------------------------------
// Line

class Line : public math::Line
{
 public:
  using math::Line::Line;
  Line(math::Line const& line) : math::Line(line) { }

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
