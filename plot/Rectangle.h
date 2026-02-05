#pragma once

#include "cairowindow/Rectangle.h"
#include "utils/Badge.h"
#include <memory>

// Forward declarations.
namespace cairowindow {

class LayerRegion;

namespace draw {
// Forward declare the draw object.
class Rectangle;
class Polyline;
} // namespace draw

template<CS> class CoordinateSystem;
template<CS> class CoordinateMapper;

} // namespace cairowindow

namespace cairowindow::plot {
// Forward declaration.
class Plot;

namespace cs {

//-----------------------------------------------------------------------------
// plot::cs::Rectangle
//
// A handle keeping a plotted Rectangle alive.
// Returned by Plot::create_rectangle(layer, rectangle_style, <args to construct a plot::cs::Rectangle>).
//
template<CS cs>
class Rectangle : public cairowindow::cs::Rectangle<cs>
{
 public:
  explicit Rectangle(cairowindow::cs::Rectangle<cs> const& rectangle) : cairowindow::cs::Rectangle<cs>(rectangle) { }
  using cairowindow::cs::Rectangle<cs>::Rectangle;

 public:
  friend class Plot;
  mutable std::shared_ptr<LayerRegion> draw_object_;       // Points to a draw::Rectangle or a draw::Polyline.

 public:
  template<typename... Args>
  void create_draw_object(utils::Badge<Plot, CoordinateSystem<cs>, CoordinateMapper<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<draw::Rectangle>(std::forward<Args>(args)...);
  }

  template<typename... Args>
  void create_polyline_draw_object(utils::Badge<Plot, CoordinateSystem<cs>, CoordinateMapper<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<draw::Polyline>(std::forward<Args>(args)...);
  }

  // Accessor for the draw object; used by Plot and CoordinateSystem.

  std::shared_ptr<LayerRegion> const& draw_object() const
  {
    return draw_object_;
  }
};

} // namespace cs

//
//-----------------------------------------------------------------------------

// The current namespace is cairowindow::plot!
//
using Rectangle = cs::Rectangle<csid::plot>;

} // namespace cairowindow::plot
