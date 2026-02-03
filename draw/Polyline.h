#pragma once

#include "LineStyle.h"
#include "cairowindow/LayerRegion.h"
#include "cairowindow/StrokeExtents.h"
#include "cairowindow/cs/Point.h"
#include <vector>
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

enum class PolylineClosedEnum {
  undefined = -1,
  open = 0,
  closed = 1
};

struct PolylineClosed {
  enum PolylineClosedEnum pc_;

  constexpr PolylineClosed(PolylineClosedEnum pc) : pc_(pc) { }
  constexpr PolylineClosed(bool closed) : pc_(closed ? PolylineClosedEnum::closed : PolylineClosedEnum::open) { }

  bool operator!=(PolylineClosedEnum pc) const { return pc_ != pc; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(pc_);
  }
#endif
};

// List the additional members of PolylineStyle.
#define cairowindow_Polyline_FOREACH_MEMBER(X, ...) \
  X(PolylineClosed, closed, PolylineClosedEnum::undefined, __VA_ARGS__)

// PolylineStyle is derived from LineStyle.
#define cairowindow_Polyline_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Line_FOREACH_STYLE_MEMBER(X, __VA_ARGS__) \
  cairowindow_Polyline_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for PolylineStyle.
struct PolylineStyleParamsDefault : LineStyleParamsDefault
{
  static constexpr PolylineClosed closed = false;
};

// Declare PolylineStyle, derived from LineStyle.
DECLARE_STYLE_WITH_BASE(Polyline, Line, PolylineStyleParamsDefault);

class Polyline : public LayerRegion
{
 private:
  std::vector<cs::Point<csid::pixels>> points_;
  PolylineStyle style_;

 public:
  Polyline(std::vector<cs::Point<csid::pixels>> points, PolylineStyle const& style) :
    points_(std::move(points)), style_(style)
  {
    ASSERT(points_.size() >= 2);
  }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::cairowindow, "draw::Polyline::do_draw(cr) [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_set_source_rgba(cr, style_.line_color().red(), style_.line_color().green(), style_.line_color().blue(), style_.line_color().alpha());
    cairo_set_line_width(cr, style_.line_width());
    if (style_.line_cap() == LineCap::round)
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
    else if (style_.line_cap() == LineCap::square)
      cairo_set_line_cap(cr, CAIRO_LINE_CAP_SQUARE);

    cairo_move_to(cr, points_[0].x(), points_[0].y());
    if (!style_.dashes().empty())
      cairo_set_dash(cr, style_.dashes().data(), style_.dashes().size(), style_.dashes_offset());
    for (std::size_t i = 1; i < points_.size(); ++i)
      cairo_line_to(cr, points_[i].x(), points_[i].y());
    if (style_.closed() != PolylineClosedEnum::open)
      cairo_close_path(cr);

    double x1, y1, x2, y2;
    cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
    cairo_stroke(cr);
    return {x1, y1, x2, y2};
  }
};

} // namespace cairowindow::draw

