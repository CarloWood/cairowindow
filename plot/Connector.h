#pragma once

#include "cairowindow/cs/Connector.h"
#include "utils/Badge.h"
#include <memory>

// Forward declarations.
namespace cairowindow {
namespace draw {
class Connector;
} // namespace draw

template<CS> class CoordinateSystem;
template<CS> class CoordinateMapper;

} // namespace cairowindow

namespace cairowindow::plot {

// Forward declaration.
class Plot;

namespace cs {

//-----------------------------------------------------------------------------
// plot::cs::Connector
//
// A handle keeping a plotted Connector alive.
// Returned by Plot::create_connector(layer, connector_style, <args to construct a cs::Connector<cs>>).
//
template<CS cs>
class Connector : public cairowindow::cs::Connector<cs>
{
 public:
  // Default constructor creates an uninitialized Connector.
  Connector() = default;
  explicit Connector(cairowindow::cs::Connector<cs> const& connector) : cairowindow::cs::Connector<cs>(connector) { }
  using cairowindow::cs::Connector<cs>::Connector;

 protected:
  mutable std::shared_ptr<draw::Connector> draw_object_;

 public:
  template<typename... Args>
  void create_draw_object(utils::Badge<Plot, cairowindow::CoordinateSystem<cs>, cairowindow::CoordinateMapper<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<draw::Connector>(std::forward<Args>(args)...);
  }

  // Erase the draw object, created with create_draw_object, if any.
  void reset()
  {
    draw_object_.reset();
  }

  // Accessor for the draw object; used by Plot and CoordinateSystem.
  std::shared_ptr<draw::Connector> const& draw_object() const
  {
    return draw_object_;
  }
};

} // namespace cs

// The current namespace is cairowindow::plot!
using Connector = cs::Connector<csid::plot>;

} // namespace cairowindow::plot
