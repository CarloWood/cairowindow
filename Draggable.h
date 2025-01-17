#pragma once

#include "utils/VectorIndex.h"
#include "debug.h"

namespace math {
class Point;
}

namespace cairowindow {
using Point = math::Point;

class Rectangle;

using ClickableIndex = utils::VectorIndex<Rectangle>;

namespace plot {

class Plot;

struct Draggable
{
  ClickableIndex index_;                // Index of this Draggable.

  virtual ~Draggable() = default;

  virtual cairowindow::Rectangle const& geometry() const = 0;
  virtual void set_position(cairowindow::Point const& new_position) = 0;
  virtual void moved(Plot* plot, cairowindow::Point const& new_position) = 0;
  virtual bool convert() const { return true; }

  void set_index(ClickableIndex index)
  {
    // Don't register the same draggable twice.
    ASSERT(index_.undefined());
    index_ = index;
  }

#ifdef CWDEBUG
  virtual void print_on(std::ostream& os) const = 0;
#endif
};

} // namespace plot
} // namespace cairowindow
