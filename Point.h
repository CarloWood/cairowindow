#pragma once

#include "Draggable.h"
#include "cs/Size.h"
#include "math/Point.h"
#include "utils/Badge.h"
#include <memory>

#include "draw/Point.h"

namespace cairowindow {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

// Forward declarations.

namespace draw {
class Point;
} // namespace draw

template<CS> class CoordinateSystem;

//--------------------------------------------------------------------------
// Current namespace is cairowindow.

namespace plot {
class Plot;

namespace cs {

// Note that adding and/or subtracting a Vector and/or Direction to/from this Point class
// returns a cairowindow::cs::Point<cs>, an object outside of namespace plot.
// The reason for that is that this object, inside namespace plot, represents a *plotted* point.
// Doing calculations with, does not magically plot the result as well.
//
template<CS cs>
class Point : public cairowindow::cs::Point<cs>, public Draggable
{
 public:
  Point() : cairowindow::cs::Point<cs>(0.0, 0.0) { }
  explicit Point(cairowindow::cs::Point<cs> const& point) : cairowindow::cs::Point<cs>(point) { }
  using cairowindow::cs::Point<cs>::Point;
  using cairowindow::cs::Point<cs>::operator=;

 protected:
  friend class Plot;
  // Drawable associated with this logical point; populated by CoordinateSystem::add_point or Plot::add_point.
  mutable std::shared_ptr<draw::Point> draw_object_;

 private:
  // Implementation of Draggable.
  cairowindow::Rectangle const& geometry() const override;
  void moved(Plot* plot, cairowindow::cs::Point<CS::plot> const& new_position) override;
  void set_position(cairowindow::cs::Point<CS::plot> const& new_position) override
  {
    this->x() = new_position.x();
    this->y() = new_position.y();
  }

 public:
  template<typename... Args>
  void create_draw_object(utils::Badge<Plot, cairowindow::CoordinateSystem<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<cairowindow::draw::Point>(std::forward<Args>(args)...);
  }

  // Accessor for the draw object; used by Plot and CoordinateSystem.

  std::shared_ptr<cairowindow::draw::Point> const& draw_object() const
  {
    return draw_object_;
  }

 public:
  void move(Plot& plot, cairowindow::cs::Point<CS::plot> const& new_position);

#ifdef CWDEBUG
 public:
  void print_on(std::ostream& os) const override { cairowindow::cs::Point<cs>::print_on(os); }
#endif
};

template<CS cs>
cairowindow::Rectangle const& Point<cs>::geometry() const
{
  return draw_object_->geometry();
}

template<CS cs>
void Point<cs>::moved(Plot* plot, cairowindow::cs::Point<CS::plot> const& new_position)
{
  // Should never call `moved` for a Point that isn't using CS::plot coordinates.
  ASSERT(false);
}

// Declare specializations of moved and move for CS::plot.

template<>
void Point<CS::plot>::moved(Plot* plot, cairowindow::cs::Point<CS::plot> const& new_position);

template<>
void Point<CS::plot>::move(Plot& UNUSED_ARG(plot), cairowindow::cs::Point<CS::plot> const& new_position);

} // namespace cs

// See remark above plot::cs::Point.
using Point = cs::Point<CS::plot>;

} // namespace plot
} // namespace cairowindow
