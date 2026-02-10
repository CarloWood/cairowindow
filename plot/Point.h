#pragma once

#include "cairowindow/Draggable.h"
#include "utils/Badge.h"
#include <memory>

// Forward declarations.
namespace cairowindow {
namespace draw {
class Point;
} // namespace draw

template<CS> class CoordinateSystem;
template<CS> class CoordinateMapper;

} // namespace cairowindow

namespace cairowindow::plot {
// Forward declaration.
class Plot;

namespace cs {

//-----------------------------------------------------------------------------
// plot::cs::Point
//
// A handle keeping a plotted Point alive.
// Returned by Plot::create_point(layer, point_style, <args to construct a plot::cs::Point<cs>>).
//
// Note that adding and/or subtracting a Vector and/or Direction to/from this Point class
// returns a math::cs::Point<cs>, an object outside of namespace plot.
// The reason for that is that this object, inside namespace plot, represents a *plotted* point.
// Doing calculations with it, does not magically plot the result as well.
//
template<CS cs>
class Point : public math::cs::Point<cs>, public Draggable
{
 public:
  // Default constructor creates a Point in the origin.
  Point() : math::cs::Point<cs>(0.0, 0.0) { }
  explicit Point(math::cs::Point<cs> const& point) : math::cs::Point<cs>(point) { }
  using math::cs::Point<cs>::Point;
  using math::cs::Point<cs>::operator=;

 protected:
  // Drawable associated with this logical point; populated by CoordinateSystem::add_point or Plot::add_point.
  mutable std::shared_ptr<draw::Point> draw_object_;

 public:
  // Implementation of Draggable.
  cairowindow::Geometry const& geometry() const override;
  void moved(math::cs::Point<csid::pixels> const& new_position_pixels) override { /*no-op*/ }

 private:
  void set_position(cairowindow::Point const& new_position) override
  {
    this->x() = new_position.x();
    this->y() = new_position.y();
  }

 public:
  template<typename... Args>
  void create_draw_object(utils::Badge<Plot, cairowindow::CoordinateSystem<cs>, cairowindow::CoordinateMapper<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<cairowindow::draw::Point>(std::forward<Args>(args)...);
  }

  // Accessor for the draw object; used by Plot and CoordinateSystem.

  std::shared_ptr<cairowindow::draw::Point> const& draw_object() const
  {
    return draw_object_;
  }

 public:
  void move_to(math::cs::Point<cs> const& new_position);

#ifdef CWDEBUG
 public:
  void print_on(std::ostream& os) const override { math::cs::Point<cs>::print_on(os); }
#endif
};

} // namespace cs

// The current namespace is cairowindow::plot!
//
// See remark above plot::cs::Point.
using Point = cs::Point<csid::plot>;

} // namespace cairowindow::plot

#include "cairowindow/draw/Point.h"
#include "cairowindow/Layer.h"

namespace cairowindow::plot::cs {

template<CS cs>
cairowindow::Geometry const& Point<cs>::geometry() const
{
  // Geometry is in csid::pixels.
  return draw_object_->geometry();
}

template<CS cs>
void Point<cs>::move_to(math::cs::Point<cs> const& new_position)
{
  // Only call move_to on points that are (registered as) draggable.
  ASSERT(!index_.undefined());
  Layer* layer = draw_object_->layer();
  Window* window = layer->window();
  window->move_draggable(this, index_, new_position);
}

} // namespace cairowindow::plot::cs
