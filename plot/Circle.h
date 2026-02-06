#pragma once

#include "cairowindow/cs/Circle.h"
#include "cairowindow/draw/Circle.h"
#include "utils/Badge.h"
#include <memory>

// Forward declarations.
namespace cairowindow {

template<CS> class CoordinateSystem;
template<CS> class CoordinateMapper;

} // namespace cairowindow

namespace cairowindow::plot {

// Forward declaration.
class Plot;

namespace cs {

//-----------------------------------------------------------------------------
// plot::cs::Circle
//
// A handle keeping a plotted Circle alive.
// Returned by Plot::create_circle(layer, circle_style, <args to construct a cs::Circle<cs>>).
//
template<CS cs>
class Circle : public cairowindow::cs::Circle<cs>
{
 public:
  Circle() = default;
  explicit Circle(cairowindow::cs::Circle<cs> const& circle) : cairowindow::cs::Circle<cs>(circle) { }
  using cairowindow::cs::Circle<cs>::Circle;

 protected:
  mutable std::shared_ptr<draw::Circle> draw_object_;

 public:
  template<typename... Args>
  void create_draw_object(utils::Badge<Plot, cairowindow::CoordinateSystem<cs>, cairowindow::CoordinateMapper<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<draw::Circle>(std::forward<Args>(args)...);
  }

  void reset()
  {
    draw_object_.reset();
  }

  std::shared_ptr<draw::Circle> const& draw_object() const
  {
    return draw_object_;
  }
};

} // namespace cs

// The current namespace is cairowindow::plot!
using Circle = cs::Circle<csid::plot>;

} // namespace cairowindow::plot

