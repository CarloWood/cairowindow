#pragma once

#include "Point.h"
#include "Vector.h"
#include "Rectangle.h"
#include "utils/VectorIndex.h"
#include "debug.h"

namespace math {
template<int N, typename T>
class Point;
}

namespace cairowindow {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

using ClickableIndex = utils::VectorIndex<Rectangle>;

namespace plot {
class Plot;

namespace cs {
template<CS cs> class Point;            // Forward declaration of plot::cs::Point.
} // namespace cs

struct Draggable
{
  ClickableIndex index_;                // Index of this Draggable.

  virtual ~Draggable() = default;

  virtual cairowindow::Geometry const& geometry() const = 0;
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
