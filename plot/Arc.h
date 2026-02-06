#pragma once

#include "cairowindow/cs/Arc.h"
#include "cairowindow/draw/Arc.h"
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
// plot::cs::Arc
//
// A handle keeping a plotted Arc alive.
// Returned by Plot::create_arc(layer, arc_style, <args to construct a cs::Arc<cs>>).
//
template<CS cs>
class Arc : public cairowindow::cs::Arc<cs>
{
 public:
  using cairowindow::cs::Arc<cs>::Arc;

  // Construct an uninitialized arc.
  Arc() = default;

 protected:
  mutable std::shared_ptr<draw::Arc> draw_object_;

 public:
  template<typename... Args>
  void create_draw_object(utils::Badge<Plot, cairowindow::CoordinateSystem<cs>, cairowindow::CoordinateMapper<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<draw::Arc>(std::forward<Args>(args)...);
  }

  void reset()
  {
    draw_object_.reset();
  }

  std::shared_ptr<draw::Arc> const& draw_object() const
  {
    return draw_object_;
  }
};

} // namespace cs

// The current namespace is cairowindow::plot!
using Arc = cs::Arc<csid::plot>;

} // namespace cairowindow::plot
