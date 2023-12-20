#pragma once

#include "Color.h"
#include "StrokeExtents.h"

namespace cairowindow {

class Layer;

class MultiRegion
{
 protected:
  Color color_;
  double line_width_;

  MultiRegion(Color const& color, double line_width) : color_(color), line_width_(line_width) { }

 public:
  virtual void draw_regions_on(Layer* layer) = 0;

 protected:
  void set_line_style(cairo_t* cr) const
  {
    cairo_set_source_rgb(cr, color_.red(), color_.green(), color_.blue());
    cairo_set_line_width(cr, line_width_);
  }

  StrokeExtents stroke(cairo_t* cr) const
  {
    DoutEntering(dc::notice, "MultiRegion::stroke(cr)");
    double x1, y1;
    double x2, y2;
    cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
    Dout(dc::notice, "cairo_stroke_extents returned " << x1 << ", " << y1 << ", " << x2 << ", " << y2);
    cairo_stroke(cr);
    return {x1, y1, x2, y2};
  }
};

} // namespace cairowindow
