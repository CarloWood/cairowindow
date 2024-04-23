#pragma once

#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow {
class Layer;

namespace draw {

class MultiRegion
{
 protected:
  Color color_;
  double line_width_;

  MultiRegion(Color const& color, double line_width) : color_(color), line_width_(line_width) { }
  virtual ~MultiRegion() = default;

 public:
  virtual void draw_regions_on(Layer* layer) = 0;

 protected:
  void set_line_style(cairo_t* cr) const
  {
    DoutEntering(dc::cairowindow, "MultiRegion::set_line_style(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_set_source_rgb(cr, color_.red(), color_.green(), color_.blue());
    cairo_set_line_width(cr, line_width_);
    cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
  }

  void set_fill_color(cairo_t* cr) const
  {
    DoutEntering(dc::cairowindow, "MultiRegion::set_fill_color(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_set_source_rgb(cr, color_.red(), color_.green(), color_.blue());
  }

  StrokeExtents stroke(cairo_t* cr) const
  {
    DoutEntering(dc::cairowindow, "MultiRegion::stroke(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    double x1, y1;
    double x2, y2;
    cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
    cairo_stroke(cr);
    return {x1, y1, x2, y2};
  }

  StrokeExtents fill(cairo_t* cr) const
  {
    DoutEntering(dc::cairowindow, "MultiRegion::fill(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    double x1, y1;
    double x2, y2;
    cairo_fill_extents(cr, &x1, &y1, &x2, &y2);
    cairo_fill(cr);
    return {x1, y1, x2, y2};
  }

#ifdef CWDEBUG
  friend std::ostream& operator<<(std::ostream& os, MultiRegion const* multi_region_ptr)
  {
    os << "MultiRegion*";
    return os;
  }
#endif
};

} // namespace draw
} // namespace cairowindow
