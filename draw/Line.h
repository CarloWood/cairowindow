#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#include "cairowindow/Style.h"
#include "cairowindow/Point.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

enum class LineCap
{
  undefined, butt, round, square
};

std::ostream& operator<<(std::ostream& os, LineCap line_cap);

#define cairowindow_Line_FOREACH_MEMBER(X, ...) \
  X(Color, line_color, Color{}, __VA_ARGS__) \
  X(double, line_width, -1.0, __VA_ARGS__) \
  X(std::vector<double>, dashes, std::vector<double>{-1.0}, __VA_ARGS__) \
  X(double, dashes_offset, 12345678.9, __VA_ARGS__) \
  X(LineCap, line_cap, LineCap::undefined, __VA_ARGS__)

#define cairowindow_Line_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Line_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for LineStyle.
struct LineStyleParamsDefault
{
  static constexpr Color line_color = color::indigo;
  static constexpr double line_width = 2.0;
  static constexpr std::vector<double> dashes = {};
  static constexpr double dashes_offset = 0.0;
  static constexpr LineCap line_cap = LineCap::butt;
};

// Declare LineStyle.
DECLARE_STYLE(Line, LineStyleParamsDefault);

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

  Line(cairowindow::Point const& point1, cairowindow::Point const& point2, LineStyle const& style) : Line(point1.x(), point1.y(), point2.x(), point2.y(), style) { }

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

} // namespace cairowindow::draw
