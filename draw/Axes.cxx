#include "sys.h"
#include "Axes.h"

namespace cairowindow::draw {

void Axes::draw_regions_on(Layer* layer)
{
  auto top_axis = [this](cairo_t* cr) -> StrokeExtents
  {
    set_line_style(cr);

    cairo_move_to(cr, geometry_.offset_x(), geometry_.offset_y());
    cairo_line_to(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y());

    return stroke(cr);
  };

  auto bottom_axis = [this](cairo_t* cr) -> StrokeExtents
  {
    set_line_style(cr);

    cairo_move_to(cr, geometry_.offset_x(), geometry_.offset_y() + geometry_.height());
    cairo_line_to(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y() + geometry_.height());

    return stroke(cr);
  };

  auto left_axis = [this](cairo_t* cr) -> StrokeExtents
  {
    set_line_style(cr);

    cairo_move_to(cr, geometry_.offset_x(), geometry_.offset_y());
    cairo_line_to(cr, geometry_.offset_x(), geometry_.offset_y() + geometry_.height());

    return stroke(cr);
  };

  auto right_axis = [this](cairo_t* cr) -> StrokeExtents
  {
    set_line_style(cr);

    cairo_move_to(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y());
    cairo_line_to(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y() + geometry_.height());

    return stroke(cr);
  };

  top_axis_.draw(layer, top_axis);
  bottom_axis_.draw(layer, bottom_axis);
  left_axis_.draw(layer, left_axis);
  right_axis_.draw(layer, right_axis);
}

} // namespace cairowindow::draw
