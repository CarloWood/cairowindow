#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

struct RectangleStyleDelta
{
  static constexpr double undefined_line_width_magic = -1.0;
  static constexpr double undefined_dashes_magic = -1.0;
  static constexpr double undefined_dashes_offset_magic = 12345678.9;

  Color line_color{};
  double line_width = undefined_line_width_magic;
  std::vector<double> dashes{undefined_dashes_magic};
  double dashes_offset = undefined_dashes_offset_magic;
};

struct RectangleStyle
{
  Color line_color = color::indigo;
  double line_width = 2.0;
  std::vector<double> dashes = {};
  double dashes_offset = 0.0;

  RectangleStyle operator()(RectangleStyleDelta delta)
  {
    RectangleStyle result{*this};
    if (delta.line_color.is_defined())
      result.line_color = delta.line_color;
    if (delta.line_width != RectangleStyleDelta::undefined_line_width_magic)
      result.line_width = delta.line_width;
    if (delta.dashes.size() != 1 || delta.dashes[0] != RectangleStyleDelta::undefined_dashes_magic)
      result.dashes = delta.dashes;
    if (delta.dashes_offset != RectangleStyleDelta::undefined_dashes_offset_magic)
      result.dashes_offset = delta.dashes_offset;
    return result;
  }
};

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
    cairo_set_source_rgba(cr, style_.line_color.red(), style_.line_color.green(), style_.line_color.blue(), style_.line_color.alpha());
    cairo_set_line_width(cr, style_.line_width);
    cairo_move_to(cr, x1_, y1_);
    if (!style_.dashes.empty())
      cairo_set_dash(cr, style_.dashes.data(), style_.dashes.size(), style_.dashes_offset);
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
