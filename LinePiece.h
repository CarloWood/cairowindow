#pragma once

#include "math/LinePiece.h"
#include <memory>

namespace cairowindow {
using LinePiece = math::LinePiece;

namespace draw {
class Line;
} // namespace draw

namespace plot {
class Plot;

//--------------------------------------------------------------------------
// LinePiece

class LinePiece : public math::LinePiece
{
 public:
  using math::LinePiece::LinePiece;
  LinePiece(math::LinePiece const& line_piece) : math::LinePiece(line_piece) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

} // namespace plot
} // namespace cairowindow
