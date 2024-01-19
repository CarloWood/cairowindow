#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

struct LineStyleDelta
{
  static constexpr double undefined_line_width_magic = -1.0;
  static constexpr double undefined_dashes_magic = -1.0;
  static constexpr double undefined_dashes_offset_magic = 12345678.9;

  Color line_color{};
  double line_width = undefined_line_width_magic;
  std::vector<double> dashes{undefined_dashes_magic};
  double dashes_offset = undefined_dashes_offset_magic;
};

struct LineStyle
{
  Color line_color = color::indigo;
  double line_width = 2.0;
  std::vector<double> dashes = {};
  double dashes_offset = 0.0;

  LineStyle operator()(LineStyleDelta delta)
  {
    LineStyle result{*this};
    if (delta.line_color.is_defined())
      result.line_color = delta.line_color;
    if (delta.line_width != LineStyleDelta::undefined_line_width_magic)
      result.line_width = delta.line_width;
    if (delta.dashes.size() != 1 || delta.dashes[0] != LineStyleDelta::undefined_dashes_magic)
      result.dashes = delta.dashes;
    if (delta.dashes_offset != LineStyleDelta::undefined_dashes_offset_magic)
      result.dashes_offset = delta.dashes_offset;
    return result;
  }
};

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

  double length() const { return std::sqrt((x2_ - x1_) * (x2_ - x1_) + (y2_ - y1_) * (y2_ - y1_)); }

 protected:
  void draw_line(cairo_t* cr)
  {
    DoutEntering(dc::cairowindow, "draw::Line(" << cr << ") [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_set_source_rgba(cr, style_.line_color.red(), style_.line_color.green(), style_.line_color.blue(), style_.line_color.alpha());
    cairo_set_line_width(cr, style_.line_width);
    cairo_move_to(cr, x1_, y1_);
    cairo_line_to(cr, x2_, y2_);
    if (!style_.dashes.empty())
      cairo_set_dash(cr, style_.dashes.data(), style_.dashes.size(), style_.dashes_offset);
    cairo_stroke(cr);
  }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::cairowindow, "draw::Line::do_draw(cr) [" << this << "]");

    draw_line(cr);

    return {std::min(x1_, x2_) - 0.5 * style_.line_width, std::min(y1_, y2_) - 0.5 * style_.line_width,
      std::max(x1_, x2_) + 0.5 * style_.line_width, std::max(y1_, y2_) + 0.5 * style_.line_width};
  }
};

} // namespace cairowindow::draw
