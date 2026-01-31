#pragma once

#include "cairowindow/Line.h"
#include "utils/Badge.h"
#include <memory>

// Forward declarations.
namespace cairowindow {

namespace draw {
// Forward declare the draw object.
class Line;
} // namespace draw

template<CS> class CoordinateSystem;

} // namespace cairowindow

namespace cairowindow::plot {
// Forward declaration.
class Plot;

namespace cs {

//-----------------------------------------------------------------------------
// plot::cs::Line
//
// A handle keeping a plotted Line alive.
// Returned by Plot::create_line(layer, line_style, <args to construct a plot::cs::Line>).
//
template<CS cs>
class Line : public cairowindow::cs::Line<cs>
{
 public:
  explicit Line(cairowindow::cs::Line<cs> const& line) : cairowindow::cs::Line<cs>(line) { }
  using cairowindow::cs::Line<cs>::Line;

  void reset()
  {
    draw_object_.reset();
  }

 protected:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;

 public:
  template<typename... Args>
  void create_draw_object(utils::Badge<Plot, cairowindow::CoordinateSystem<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<draw::Line>(std::forward<Args>(args)...);
  }

  // Accessor for the draw object; used by Plot and CoordinateSystem.

  std::shared_ptr<draw::Line> const& draw_object() const
  {
    return draw_object_;
  }
};

} // namespace cs

//
//-----------------------------------------------------------------------------

// The current namespace is cairowindow::plot!
//
using Line = cs::Line<CS::plot>;

} // namespace cairowindow::plot
