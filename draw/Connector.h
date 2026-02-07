#pragma once

#include "Line.h"
#include "Shape.h"
#include "ArrowHead.h"
#include "cairowindow/cs/Connector.h"
#include "cairowindow/Layer.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#endif

namespace cairowindow::draw {

// List the additional members of ConnectorStyle.
#define cairowindow_Connector_FOREACH_MEMBER(X, ...) \
  X(ArrowHeadStyle, arrow_head_from, ArrowHeadStyle{{.position = undefined_shape_position}}, __VA_ARGS__) \
  X(ArrowHeadStyle, arrow_head_to, ArrowHeadStyle{{.position = undefined_shape_position}}, __VA_ARGS__)

// ConnectorStyle is derived from LineStyle.
#define cairowindow_Connector_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Line_FOREACH_STYLE_MEMBER(X, __VA_ARGS__) \
  cairowindow_Connector_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for ConnectorStyle.
struct ConnectorStyleParamsDefault : LineStyleParamsDefault
{
  static ArrowHeadStyle const arrow_head_from;
  static ArrowHeadStyle const arrow_head_to;
};

inline bool operator!=(ArrowHeadStyle const& style1, ArrowHeadStyle const& style2)
{
  // Only use this operator to compare with the undefined magic.
  ASSERT(style2.position() == undefined_shape_position);
  return style1.position() != undefined_shape_position;
}

// Declare ConnectorStyle, derived from LineStyle.
DECLARE_STYLE_WITH_BASE(Connector, Line, ConnectorStyleParamsDefault);

class Connector : public Line
{
 public:
  using Direction = math::cs::Direction<csid::pixels>;
  using ArrowHeadShape = cairowindow::cs::Connector<csid::plot>::ArrowHeadShape;
  static ShapeEnum to_ShapeEnum(ArrowHeadShape arrow_head_shape) { return static_cast<ShapeEnum>(arrow_head_shape + number_of_shapes); }

 private:
  std::shared_ptr<ArrowHead> arrow_head_from_;
  std::shared_ptr<ArrowHead> arrow_head_to_;

 public:
  Connector(double x1, double y1, double x2, double y2, ConnectorStyle const& connector_style,
      ArrowHeadShape arrow_header_shape_from, ArrowHeadShape arrow_header_shape_to) :
    Line(x1, y1, x2, y2, connector_style)
  {
    // Don't draw arrow heads when the length of the Connector is near zero.
    if (std::abs(x2 - x1) < 0.01 && std::abs(y2 - y1) < 0.01)
      return;
    if (arrow_header_shape_from != cairowindow::cs::Connector<csid::plot>::no_arrow)
    {
      Color line_color = connector_style.arrow_head_from().line_color();
      if (line_color.is_transparent())
        line_color = connector_style.line_color();
      Color fill_color = connector_style.arrow_head_from().fill_color();
      if (fill_color.is_transparent())
        fill_color = connector_style.line_color();
      arrow_head_from_ =
        std::make_shared<ArrowHead>(x1, y1, Direction{{x2, y2}, {x1, y1}}, // Pointing to the tip at x1,y1
            connector_style.arrow_head_from()({.line_color = line_color, .fill_color = fill_color,
              .shape = to_ShapeEnum(arrow_header_shape_from)}));
    }
    if (arrow_header_shape_to != cairowindow::cs::Connector<csid::plot>::no_arrow)
    {
      Color line_color = connector_style.arrow_head_to().line_color();
      if (line_color.is_transparent())
        line_color = connector_style.line_color();
      Color fill_color = connector_style.arrow_head_to().fill_color();
      if (fill_color.is_transparent())
        fill_color = connector_style.line_color();
      arrow_head_to_ = std::make_shared<ArrowHead>(x2, y2, Direction{{x1, y1}, {x2, y2}},
          connector_style.arrow_head_to()({.line_color = line_color, .fill_color = fill_color,
            .shape = to_ShapeEnum(arrow_header_shape_to)}));
    }
  }

  void draw_arrow_heads(boost::intrusive_ptr<Layer> const& layer);

  void update_from(double overshoot);
  void update_to(double overshoot);

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::cairowindow, "draw::Connector::do_draw(cr) [" << this << "]");

    draw_line(cr);

    return {std::min(x1_, x2_) - 0.5 * style_.line_width(), std::min(y1_, y2_) - 0.5 * style_.line_width(),
      std::max(x1_, x2_) + 0.5 * style_.line_width(), std::max(y1_, y2_) + 0.5 * style_.line_width()};
  }
};

} // namespace cairowindow::draw
