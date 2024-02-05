#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

// List the additional members of RectangleStyle (none).
#define cairowindow_Rectangle_FOREACH_MEMBER(X, ...)

// RectangleStyle is derived from LineStyle.
#define cairowindow_Rectangle_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Line_FOREACH_MEMBER(X, __VA_ARGS__) \
  cairowindow_Rectangle_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for RectangleStyle.
struct RectangleStyleParamsDefault : LineStyleParamsDefault
{
};

// Declare RectangleStyle, derived from LineStyle.
DECLARE_STYLE_WITH_BASE(Rectangle, Line, RectangleStyleParamsDefault);

class Rectangle : public LayerRegion
{
 private:
  double x1_;
  double y1_;
  double x2_;
  double y2_;
  RectangleStyle style_;

 public:
  Rectangle(double x1, double y1, double x2, double y2, RectangleStyle const& style) :
    x1_(x1), y1_(y1), x2_(x2), y2_(y2), style_(style) { }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::cairowindow, "draw::Rectangle::do_draw(cr) [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_set_source_rgba(cr, style_.line_color().red(), style_.line_color().green(), style_.line_color().blue(), style_.line_color().alpha());
    cairo_set_line_width(cr, style_.line_width());
    cairo_move_to(cr, x1_, y1_);
    if (!style_.dashes().empty())
      cairo_set_dash(cr, style_.dashes().data(), style_.dashes().size(), style_.dashes_offset());
    cairo_line_to(cr, x2_, y1_);
    cairo_line_to(cr, x2_, y2_);
    cairo_line_to(cr, x1_, y2_);
    cairo_close_path(cr);
    double x1, y1, x2, y2;
    cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
    cairo_stroke(cr);

    return {x1, y1, x2, y2};
  }
};

} // namespace cairowindow::draw
