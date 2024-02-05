#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

// List the additional members of ArcStyle (none).
#define cairowindow_Arc_FOREACH_MEMBER(X, ...)

// ArcStyle is derived from LineStyle.
#define cairowindow_Arc_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Line_FOREACH_MEMBER(X, __VA_ARGS__) \
  cairowindow_Arc_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for ArcStyle.
struct ArcStyleParamsDefault : LineStyleParamsDefault
{
};

// Declare ArcStyle, derived from LineStyle.
DECLARE_STYLE_WITH_BASE(Arc, Line, ArcStyleParamsDefault);

class Arc : public LayerRegion
{
 protected:
  double center_x_;
  double center_y_;
  double pixel_start_angle_;
  double pixel_end_angle_;
  double radius_;
  ArcStyle style_;

 public:
  Arc(double center_x, double center_y, double pixel_start_angle, double pixel_end_angle, double radius, ArcStyle const& style) :
    center_x_(center_x), center_y_(center_y), pixel_start_angle_(pixel_start_angle), pixel_end_angle_(pixel_end_angle), radius_(radius), style_(style) { }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::cairowindow, "draw::Arc::do_draw(cr) [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_set_source_rgba(cr, style_.line_color().red(), style_.line_color().green(), style_.line_color().blue(), style_.line_color().alpha());
    cairo_set_line_width(cr, style_.line_width());
    if (!style_.dashes().empty())
      cairo_set_dash(cr, style_.dashes().data(), style_.dashes().size(), style_.dashes_offset());
    cairo_arc(cr, center_x_, center_y_, radius_, pixel_start_angle_, pixel_end_angle_);
    double x1, y1, x2, y2;
    cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
    cairo_stroke(cr);

    return {x1, y1, x2, y2};
  }
};

} // namespace cairowindow::draw
