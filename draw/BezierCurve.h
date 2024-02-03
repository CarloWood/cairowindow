#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

struct BezierCurveStyle
{
  Color line_color = color::gray;
  double line_width = 2.0;
  std::vector<double> dashes = {};
  double dashes_offset = 0.0;
};

class BezierCurve : public LayerRegion
{
 protected:
  double x0_;
  double y0_;
  double x1_;
  double y1_;
  double x2_;
  double y2_;
  double x3_;
  double y3_;
  BezierCurveStyle style_;

 public:
  BezierCurve(double x0, double y0, double x1, double y1, double x2, double y2, double x3, double y3,
      BezierCurveStyle const& bezier_curve_style) :
    x0_(x0), y0_(y0), x1_(x1), y1_(y1), x2_(x2), y2_(y2), x3_(x3), y3_(y3), style_(bezier_curve_style) { }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::cairowindow, "draw::BezierCurve::do_draw(cr) [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_set_source_rgba(cr, style_.line_color.red(), style_.line_color.green(), style_.line_color.blue(), style_.line_color.alpha());
    cairo_set_line_width(cr, style_.line_width);
    if (!style_.dashes.empty())
      cairo_set_dash(cr, style_.dashes.data(), style_.dashes.size(), style_.dashes_offset);

    cairo_move_to(cr, x0_, y0_);
    cairo_curve_to(cr, x1_, y1_, x2_, y2_, x3_, y3_);

    double x1, y1, x2, y2;
    cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
    cairo_stroke(cr);

    return {x1, y1, x2, y2};
  }
};

} // namespace cairowindow::draw
