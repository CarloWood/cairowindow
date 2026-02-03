#ifndef CAIROWINDOW_DRAW_LINE_H
#define CAIROWINDOW_DRAW_LINE_H

#include "LineStyle.h"
#include "cairowindow/LayerRegion.h"
#include "cairowindow/StrokeExtents.h"
#include "cairowindow/cs/Point.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow {
using Point = cs::Point<csid::plot>;

namespace draw {

class Line : public LayerRegion
{
 protected:
  double x1_;
  double y1_;
  double x2_;
  double y2_;
  LineStyle style_;

 public:
  Line(double x1, double y1, double x2, double y2, LineStyle const& style) : x1_(x1), y1_(y1), x2_(x2), y2_(y2), style_(style)
  {
    ASSERT(!std::isnan(x1) && !std::isnan(y1) && !std::isnan(x2) && !std::isnan(y2));
  }

  inline Line(cairowindow::cs::Point<csid::pixels> const& point1, cairowindow::cs::Point<csid::pixels> const& point2, LineStyle const& style);

  double length() const { return std::sqrt((x2_ - x1_) * (x2_ - x1_) + (y2_ - y1_) * (y2_ - y1_)); }

 protected:
  void draw_line(cairo_t* cr)
  {
    DoutEntering(dc::cairowindow, "draw::Line(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_set_source_rgba(cr, style_.line_color().red(), style_.line_color().green(), style_.line_color().blue(), style_.line_color().alpha());
    cairo_set_line_width(cr, style_.line_width());
    if (style_.line_cap() == LineCap::round)
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
    else if (style_.line_cap() == LineCap::square)
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);
    cairo_move_to(cr, x1_, y1_);
    cairo_line_to(cr, x2_, y2_);
    if (!style_.dashes().empty())
      cairo_set_dash(cr, style_.dashes().data(), style_.dashes().size(), style_.dashes_offset());
    cairo_stroke(cr);
  }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::cairowindow, "draw::Line::do_draw(cr) [" << this << "]");

    draw_line(cr);

    return {std::min(x1_, x2_) - 0.5 * style_.line_width(), std::min(y1_, y2_) - 0.5 * style_.line_width(),
      std::max(x1_, x2_) + 0.5 * style_.line_width(), std::max(y1_, y2_) + 0.5 * style_.line_width()};
  }
};

} // namespace draw
} // namespace cairowindow

#endif // CAIROWINDOW_DRAW_LINE_H

#include "Point.h"

#ifndef CAIROWINDOW_DRAW_LINE_H_definitions
#define CAIROWINDOW_DRAW_LINE_H_definitions

namespace cairowindow::draw {

//inline
Line::Line(cairowindow::cs::Point<csid::pixels> const& point1, cairowindow::cs::Point<csid::pixels> const& point2, LineStyle const& style) :
  Line(point1.x(), point1.y(), point2.x(), point2.y(), style)
{
}

} // namespace cairowindow::draw

#endif // CAIROWINDOW_DRAW_LINE_H_definitions
