#pragma once

#include "cairowindow/LinePiece.h"
#include <memory>

namespace cairowindow::draw {
// Forward declare the draw object.
class Line;
} // namespace cairowindow::draw

namespace cairowindow::plot {
// Forward declaration.
class Plot;

//-----------------------------------------------------------------------------
// plot::LinePiece
//
// A handle keeping a plotted LinePiece alive.
// Returned by Plot::create_line(layer, line_style, [line_extend,] <args to construct a plot::LinePiece>).
//
class LinePiece : public cairowindow::LinePiece
{
 public:
  explicit LinePiece(cairowindow::LinePiece const& line_piece) : cairowindow::LinePiece(line_piece) { }
  using cairowindow::LinePiece::LinePiece;

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

//
//-----------------------------------------------------------------------------

} // namespace cairowindow::plot
