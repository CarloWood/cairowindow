#pragma once

#include "cs/Line.h"
#include <memory>

namespace cairowindow {
using Line = cs::Line<CS::plot>;

namespace draw {
class Line;
} // namespace draw

namespace plot {
class Plot;

//--------------------------------------------------------------------------
// Line

class Line : public cairowindow::cs::Line<CS::plot>
{
 public:
  explicit Line(cairowindow::cs::Line<CS::plot> const& line) : cairowindow::cs::Line<CS::plot>(line) { }
  using cairowindow::cs::Line<CS::plot>::Line;

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
