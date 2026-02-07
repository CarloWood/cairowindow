#pragma once

#include "cairowindow/draw/Line.h"
#include "math/cs/LinePiece.h"
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
// plot::cs::LinePiece
//
// A handle keeping a plotted LinePiece alive.
// Returned by Plot::create_line(layer, line_style, [line_extend,] <args to construct a cs::LinePiece<cs>>).
//
template<CS cs>
class LinePiece : public math::cs::LinePiece<cs>
{
 public:
  // Default constructor creates an uninitialized LinePiece.
  LinePiece() = default;
  explicit LinePiece(math::cs::LinePiece<cs> const& line_piece) : math::cs::LinePiece<cs>(line_piece) { }
  using math::cs::LinePiece<cs>::LinePiece;

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
using LinePiece = cs::LinePiece<csid::plot>;

} // namespace cairowindow::plot
