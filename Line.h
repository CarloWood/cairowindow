#pragma once

#include "Direction.h"
#include "Point.h"
#include <memory>

namespace cairowindow {

class Line
{
 private:
  Point point_;
  Direction direction_;

 public:
  // Construct an undefined line.
  Line() = default;
  // Construct a line through point with direction.
  Line(Point const& point, Direction const& direction) : point_(point), direction_(direction) { }

  Point const& point() const { return point_; }
  Direction const& direction() const { return direction_; }

  operator Direction const&() const { return direction_; }

  Point intersection_with(Line const& line2) const;
};

namespace draw {
class Line;
} // namespace draw

namespace plot {
class Plot;

//--------------------------------------------------------------------------
// Line

class Line : public cairowindow::Line
{
 public:
  using cairowindow::Line::Line;
  Line(cairowindow::Line const& line) : cairowindow::Line(line) { }

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
