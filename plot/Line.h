#pragma once

#include "cairowindow/Line.h"
#include "utils/Badge.h"
#include <memory>

// Forward declarations.
namespace cairowindow {
namespace draw {
class Line;
} // namespace draw

template<CS> class CoordinateSystem;
template<CS> class CoordinateMapper;

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
  // Default constructor creates an uninitialized Line.
  Line() = default;
  explicit Line(cairowindow::cs::Line<cs> const& line) : cairowindow::cs::Line<cs>(line) { }
  using cairowindow::cs::Line<cs>::Line;

 protected:
  mutable std::shared_ptr<draw::Line> draw_object_;

 public:
  template<typename... Args>
  void create_draw_object(utils::Badge<Plot, cairowindow::CoordinateSystem<cs>, cairowindow::CoordinateMapper<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<draw::Line>(std::forward<Args>(args)...);
  }

  // Erase the draw object, created with create_draw_object, if any.
  void reset()
  {
    draw_object_.reset();
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
using Line = cs::Line<csid::plot>;

} // namespace cairowindow::plot
