#pragma once

#include "Draggable.h"
#include "math/Point.h"
#include <memory>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace cairowindow {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

using Point = cairowindow::Point;

namespace draw {
class Point;
} // namespace draw

namespace plot {
class Plot;

class Point : public cairowindow::Point, public Draggable
{
 public:
  Point() = default;
  Point(cairowindow::Point const& point) : cairowindow::Point(point) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Point> draw_object_;

  // Implementation of Draggable.
  cairowindow::Rectangle const& geometry() const override;
  void moved(Plot* plot, cairowindow::Point const& new_position) override;
  void set_position(cairowindow::Point const& new_position) override
  {
    x_ = new_position.x();
    y_ = new_position.y();
  }

 public:
  void move(Plot& plot, cairowindow::Point const& new_position);

#ifdef CWDEBUG
 public:
  void print_on(std::ostream& os) const override { cairowindow::Point::print_on(os); }
#endif
};

} // namespace plot
} // namespace cairowindow
