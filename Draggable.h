#pragma once

#include "cs/Point.h"
#include "cs/Vector.h"
#include "utils/VectorIndex.h"
#include "debug.h"

namespace math {
template<int N, typename T>
class Point;
}

namespace cairowindow {
// Defined here - because this header is included from Plot.h; the most basis types that a user works with in that case are in 'plot' coordinates.
using Point = cs::Point<CS::plot>;
using Vector = cs::Vector<CS::plot>;
using Direction = cs::Direction<CS::plot>;

class Rectangle;

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
