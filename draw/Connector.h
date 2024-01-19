#pragma once

#include "Line.h"
#include "Shape.h"
#include "ArrowHead.h"
#include "cairowindow/Connector.h"
#include "cairowindow/Layer.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#endif

namespace cairowindow::draw {

class Connector : public Line
{
  using ArrowHeadShape = cairowindow::Connector::ArrowHeadShape;
  static ShapeEnum shape(ArrowHeadShape arrow_head_shape) { return static_cast<ShapeEnum>(arrow_head_shape + number_of_shapes); }

 private:
  ArrowHeadShape arrow_head_shape_from_;
  ArrowHeadShape arrow_head_shape_to_;
  std::shared_ptr<ArrowHead> arrow_head_from_;
  std::shared_ptr<ArrowHead> arrow_head_to_;

 public:
  Connector(double x1, double y1, double x2, double y2, LineStyle const& style, Color fill_color,
      cairowindow::Connector::ArrowHeadShape arrow_head_shape_from, cairowindow::Connector::ArrowHeadShape arrow_head_shape_to) :
    Line(x1, y1, x2, y2, style), arrow_head_shape_from_(arrow_head_shape_from), arrow_head_shape_to_(arrow_head_shape_to)
  {
    // Don't draw arrow heads when the length of the Connector is near zero.
    if (std::abs(x2 - x1) < 0.01 && std::abs(y2 - y1) < 0.01)
    {
      arrow_head_shape_from_ = cairowindow::Connector::no_arrow;
      arrow_head_shape_to_ = cairowindow::Connector::no_arrow;
      return;
    }
    if (arrow_head_shape_from != cairowindow::Connector::no_arrow)
      arrow_head_from_ = std::make_shared<ArrowHead>(x1, y1, Direction{{x2, y2}, {x1, y1}}, // Pointing to the tip at x1,y1
          ArrowHeadStyle{.line_color = style.line_color, .fill_color = fill_color, .line_width = style.line_width,
           .shape = Connector::shape(arrow_head_shape_from)});
    if (arrow_head_shape_to != cairowindow::Connector::no_arrow)
      arrow_head_to_ = std::make_shared<ArrowHead>(x2, y2, Direction{{x1, y1}, {x2, y2}},
          ArrowHeadStyle{.line_color = style.line_color, .fill_color = fill_color, .line_width = style.line_width,
           .shape = Connector::shape(arrow_head_shape_to)});
  }

  ArrowHeadShape arrow_head_shape_from() const { return arrow_head_shape_from_; }
  ArrowHeadShape arrow_head_shape_to() const { return arrow_head_shape_to_; }

  void draw_arrow_heads(boost::intrusive_ptr<Layer> const& layer);

  void update_from(double overshoot);
  void update_to(double overshoot);

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::cairowindow, "draw::Connector::do_draw(cr) [" << this << "]");

    draw_line(cr);

    return {std::min(x1_, x2_) - 0.5 * style_.line_width, std::min(y1_, y2_) - 0.5 * style_.line_width,
      std::max(x1_, x2_) + 0.5 * style_.line_width, std::max(y1_, y2_) + 0.5 * style_.line_width};
  }
};

} // namespace cairowindow::draw
