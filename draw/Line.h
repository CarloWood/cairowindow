#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"

namespace cairowindow::draw {

struct LineStyle
{
  Color line_color = color::indigo;
  double line_width = 2.0;
  std::vector<double> dashes = {};
  double dashes_offset = 0.0;
};

class Line : public LayerRegion
{
 private:
  double x1_;
  double y1_;
  double x2_;
  double y2_;
  LineStyle style_;

 public:
  Line(double x1, double y1, double x2, double y2, LineStyle style) : x1_(x1), y1_(y1), x2_(x2), y2_(y2), style_(style) { }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::notice, "draw::Line::do_draw(cr) [" << this << "]");

    cairo_set_source_rgba(cr, style_.line_color.red(), style_.line_color.green(), style_.line_color.blue(), style_.line_color.alpha());
    cairo_set_line_width(cr, style_.line_width);
    cairo_move_to(cr, x1_, y1_);
    cairo_line_to(cr, x2_, y2_);
    if (!style_.dashes.empty())
      cairo_set_dash(cr, style_.dashes.data(), style_.dashes.size(), style_.dashes_offset);
    cairo_stroke(cr);
    return {std::min(x1_, x2_) - 0.5 * style_.line_width, std::min(y1_, y2_) - 0.5 * style_.line_width,
      std::max(x1_, x2_) + 0.5 * style_.line_width, std::max(y1_, y2_) + 0.5 * style_.line_width};
  }
};

} // namespace cairowindow::draw
