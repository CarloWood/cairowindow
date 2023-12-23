#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"

namespace cairowindow::draw {

struct LineStyle
{
  Color line_color = color::indigo;
  double line_width = 2.0;
};

class Line : public LayerRegion
{
  // Lets not get confused with draw::Rectangle (in case that is #include-d).
  using Rectangle = cairowindow::Rectangle;

 private:
  Rectangle geometry_;
  LineStyle style_;

 public:
  Line(Rectangle const& geometry, LineStyle style) : geometry_(geometry), style_(style) { }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::notice, "draw::Line::do_draw(cr) [" << this << "]");

    cairo_set_source_rgba(cr, style_.line_color.red(), style_.line_color.green(), style_.line_color.blue(), style_.line_color.alpha());
    cairo_set_line_width(cr, style_.line_width);
    cairo_move_to(cr, geometry_.offset_x(),                     geometry_.offset_y());
    cairo_line_to(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y() + geometry_.height());
    cairo_stroke(cr);
    return {geometry_.offset_x() - 0.5 * style_.line_width, geometry_.offset_y() - 0.5 * style_.line_width,
      geometry_.offset_x() + geometry_.width() + 0.5 * style_.line_width, geometry_.offset_y() + geometry_.height() + 0.5 * style_.line_width};
  }
};

} // namespace cairowindow::draw
