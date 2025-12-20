#pragma once

#include "LinePiece.h"
#include "draw/Shape.h"
#include "utils/square.h"
#include <cmath>
#include <memory>

namespace cairowindow {

class Connector : public LinePiece
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
  ArrowHeadShape arrow_head_shape_from_;
  ArrowHeadShape arrow_head_shape_to_;

 public:
  Connector() = default;

  Connector(Point const& from, Point const& to, ArrowHeadShape arrow_head_shape_from, ArrowHeadShape arrow_head_shape_to) :
    LinePiece(from, to), arrow_head_shape_from_(arrow_head_shape_from), arrow_head_shape_to_(arrow_head_shape_to) { }

  Connector(Point const& from, Point const& to, ArrowHeadShape arrow_head_shape_to = open_arrow) :
    LinePiece(from, to), arrow_head_shape_from_(no_arrow), arrow_head_shape_to_(arrow_head_shape_to) { }

  ArrowHeadShape arrow_head_shape_from() const { return arrow_head_shape_from_; }
  ArrowHeadShape arrow_head_shape_to() const { return arrow_head_shape_to_; }
};

namespace draw {
class Connector;
} // namespace draw

namespace plot {
class Plot;

class Connector : public cairowindow::Connector
{
 public:
  explicit Connector(cairowindow::Connector const& connector) : cairowindow::Connector(connector) { }
  using cairowindow::Connector::Connector;

  void reset()
  {
    draw_object_.reset();
  }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Connector> draw_object_;
};

} // namespace plot
} // namespace cairowindow
