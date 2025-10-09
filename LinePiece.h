#pragma once

#include "math/LinePiece.h"
#include <memory>

namespace cairowindow {
using LinePiece = math::LinePiece<2>;

namespace draw {
class Line;
} // namespace draw

namespace plot {
class Plot;

//--------------------------------------------------------------------------
// LinePiece

class LinePiece : public math::LinePiece<2>
{
 public:
  using math::LinePiece<2>::LinePiece;
  LinePiece(math::LinePiece<2> const& line_piece) : math::LinePiece<2>(line_piece) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

} // namespace plot
} // namespace cairowindow
