#pragma once

#include "LayerRegion.h"
#include "Color.h"

namespace cairowindow::draw {

class Line : public LayerRegion
{
 private:
  Rectangle geometry_;
  Color color_;
  double line_width_;

 public:
  Line(Layer* layer, Rectangle const& geometry, Color const& color, double line_width = 2.0) :
    LayerRegion(layer), geometry_(geometry), color_(color), line_width_(line_width) { }

 private:
  Rectangle do_draw(cairo_t* cr) override
  {
    cairo_set_source_rgb(cr, color_.red(), color_.green(), color_.blue());
    cairo_set_line_width(cr, line_width_);
    cairo_move_to(cr, geometry_.offset_x(),                     geometry_.offset_y());
    cairo_line_to(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y() + geometry_.height());
    cairo_stroke(cr);
    return {geometry_.offset_x() - 0.5 * line_width_, geometry_.offset_y() - 0.5 * line_width_,
            geometry_.width() + line_width_, geometry_.height() + line_width_};
  }
};

} // namespace cairowindow::draw
