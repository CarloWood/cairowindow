#pragma once

#include "utils/VectorIndex.h"

namespace cairowindow {

class Rectangle;
class Point;

using ClickableIndex = utils::VectorIndex<Rectangle>;

namespace plot {

class Plot;

struct Draggable
{
  ClickableIndex index_;                // Index of this Draggable.

  virtual ~Draggable() = default;

  virtual Rectangle const& geometry() const = 0;
  virtual void moved(Plot* plot, cairowindow::Point const& new_position) = 0;
  virtual bool convert() const { return true; }

  void set_index(ClickableIndex index) { index_ = index; }

#ifdef CWDEBUG
  virtual void print_on(std::ostream& os) const = 0;
#endif
};

} // namespace plot
} // namespace cairowindow
