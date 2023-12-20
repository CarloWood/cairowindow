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

  region_.set_layer(layer);
  region_.draw(top_axis);
}

} // namespace cairowindow::draw
