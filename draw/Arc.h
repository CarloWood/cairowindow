#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"
#include "cairowindow/StrokeExtents.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

struct ArcStyleDelta
{
  static constexpr double undefined_line_width_magic = -1.0;
  static constexpr double undefined_dashes_magic = -1.0;
  static constexpr double undefined_dashes_offset_magic = 12345678.9;

  Color line_color{};
  double line_width = undefined_line_width_magic;
  std::vector<double> dashes{undefined_dashes_magic};
  double dashes_offset = undefined_dashes_offset_magic;
};

struct ArcStyle
{
  Color line_color = color::indigo;
  double line_width = 2.0;
  std::vector<double> dashes = {};
  double dashes_offset = 0.0;

  ArcStyle operator()(ArcStyleDelta delta)
  {
    ArcStyle result{*this};
    if (delta.line_color.is_defined())
      result.line_color = delta.line_color;
    if (delta.line_width != ArcStyleDelta::undefined_line_width_magic)
      result.line_width = delta.line_width;
    if (delta.dashes.size() != 1 || delta.dashes[0] != ArcStyleDelta::undefined_dashes_magic)
      result.dashes = delta.dashes;
    if (delta.dashes_offset != ArcStyleDelta::undefined_dashes_offset_magic)
      result.dashes_offset = delta.dashes_offset;
    return result;
  }
};

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
    cairo_set_source_rgba(cr, style_.line_color.red(), style_.line_color.green(), style_.line_color.blue(), style_.line_color.alpha());
    cairo_set_line_width(cr, style_.line_width);
    cairo_arc(cr, center_x_, center_y_, radius_, pixel_start_angle_, pixel_end_angle_);
    if (!style_.dashes.empty())
      cairo_set_dash(cr, style_.dashes.data(), style_.dashes.size(), style_.dashes_offset);
    double x1, y1, x2, y2;
    cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
    cairo_stroke(cr);

    return {x1, y1, x2, y2};
  }
};

} // namespace cairowindow::draw
