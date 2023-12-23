#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"

namespace cairowindow::draw {

struct RectangleStyle
{
  Color line_color = color::red;
  Color fill_color = color::transparent;
  double line_width = 2.0;
};

class Rectangle : public LayerRegion
{
 private:
  cairowindow::Rectangle geometry_;
  RectangleStyle style_;

 public:
  Rectangle(cairowindow::Rectangle const& geometry, RectangleStyle style) : geometry_(geometry), style_(style) { }

 private:
  cairowindow::StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::notice, "draw::Rectangle::do_draw(cr) [" << this << "]");

    bool do_stroke = !style_.line_color.is_transparent();
    bool do_fill = !style_.fill_color.is_transparent();
    ASSERT(do_stroke || do_fill);

    cairo_move_to(cr, geometry_.offset_x(),                     geometry_.offset_y());
    cairo_line_to(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y());
    cairo_line_to(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y() + geometry_.height());
    cairo_line_to(cr, geometry_.offset_x(),                     geometry_.offset_y() + geometry_.height());
    cairo_line_to(cr, geometry_.offset_x(),                     geometry_.offset_y());

    double x1, y1;
    double x2, y2;

    if (do_fill)
      cairo_set_source_rgba(cr, style_.fill_color.red(), style_.fill_color.green(), style_.fill_color.blue(), style_.fill_color.alpha());
    if (do_stroke)
    {
      if (do_fill)
        cairo_fill_preserve(cr);
      cairo_set_source_rgba(cr, style_.line_color.red(), style_.line_color.green(), style_.line_color.blue(), style_.line_color.alpha());
      cairo_set_line_width(cr, style_.line_width);
      cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
      cairo_stroke(cr);
    }
    else
    {
      cairo_fill_extents(cr, &x1, &y1, &x2, &y2);
      cairo_fill(cr);
    }

    return {x1, y1, x2, y2};
  }
};

} // namespace cairowindow::draw
