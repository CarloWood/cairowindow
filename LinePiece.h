#pragma once

#include "cs/LinePiece.h"
#include <memory>

namespace cairowindow {
using LinePiece = cs::LinePiece<CS::plot>;

namespace draw {
class Line;
} // namespace draw

namespace plot {
class Plot;

//--------------------------------------------------------------------------
// LinePiece

class LinePiece : public cairowindow::LinePiece
{
 public:
  explicit LinePiece(cairowindow::LinePiece const& line_piece) : cairowindow::LinePiece(line_piece) { }
  using cairowindow::LinePiece::LinePiece;

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

} // namespace plot
} // namespace cairowindow
