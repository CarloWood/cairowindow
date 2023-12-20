#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"

namespace cairowindow::draw {

class Line : public LayerRegion
{
  // Lets not get confused with draw::Rectangle (in case that is #include-d).
  using Rectangle = cairowindow::Rectangle;

 private:
  Rectangle geometry_;
  Color color_;
  double line_width_;

 public:
  Line(Rectangle const& geometry, Color const& color, double line_width = 2.0) :
    geometry_(geometry), color_(color), line_width_(line_width) { }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::notice, "draw::Line::do_draw(cr) [" << this << "]");

    cairo_set_source_rgb(cr, color_.red(), color_.green(), color_.blue());
    cairo_set_line_width(cr, line_width_);
    cairo_move_to(cr, geometry_.offset_x(),                     geometry_.offset_y());
    cairo_line_to(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y() + geometry_.height());
    cairo_stroke(cr);
    return {geometry_.offset_x() - 0.5 * line_width_, geometry_.offset_y() - 0.5 * line_width_,
      geometry_.offset_x() + geometry_.width() + 0.5 * line_width_, geometry_.offset_y() + geometry_.height() + 0.5 * line_width_};
  }
};

} // namespace cairowindow::draw
