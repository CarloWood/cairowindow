#pragma once

#include "math/cs/LinePiece.h"

namespace cairowindow::cs {

template<CS cs>
class Connector : public LinePiece<cs>
{
 public:
  enum ArrowHeadShape
  {
    no_arrow,
    open_arrow,
    filled_arrow,
    diamond_arrow,
    circle_arrow
  };

 private:
  ArrowHeadShape arrow_head_shape_from_{no_arrow};
  ArrowHeadShape arrow_head_shape_to_{open_arrow};

 public:
  Connector() = default;

  Connector(Point<cs> const& from, Point<cs> const& to, ArrowHeadShape arrow_head_shape_from, ArrowHeadShape arrow_head_shape_to) :
    LinePiece<cs>(from, to), arrow_head_shape_from_(arrow_head_shape_from), arrow_head_shape_to_(arrow_head_shape_to) { }

  Connector(Point<cs> const& from, Point<cs> const& to, ArrowHeadShape arrow_head_shape_to = open_arrow) :
    LinePiece<cs>(from, to), arrow_head_shape_from_(no_arrow), arrow_head_shape_to_(arrow_head_shape_to) { }

  ArrowHeadShape arrow_head_shape_from() const { return arrow_head_shape_from_; }
  ArrowHeadShape arrow_head_shape_to() const { return arrow_head_shape_to_; }
};

} // namespace cairowindow::cs
