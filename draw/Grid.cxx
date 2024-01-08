#include "sys.h"
#include "Grid.h"
#include "PlotArea.h"
#include "cairowindow/Layer.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#endif

namespace cairowindow::draw {

void Grid::draw_regions_on(Layer* layer)
{
  for (int axis = plot::x_axis; axis <= plot::y_axis; ++axis)
  {
    double x1 = geometry_.offset_x();
    double y1 = geometry_.offset_y() + geometry_.height();
    double delta = ((axis == plot::x_axis) ? geometry_.width() : geometry_.height()) / ticks_[axis];
    double x2 = x1;
    double y2 = y1;
    if (axis == plot::x_axis)
      y2 -= geometry_.height();
    else
      x2 += geometry_.width();
    Dout(dc::cairowindow, "axis = " << axis << "; ticks_[axis] = " << ticks_[axis]);
    for (int tick = 1; tick < ticks_[axis]; ++tick)
    {
      if (axis == plot::x_axis)
      {
        x1 += delta;
        x2 += delta;
      }
      else
      {
        y1 -= delta;
        y2 -= delta;
      }
      lines_.emplace_back(std::make_shared<Line>(x1, y2, x2, y1, LineStyle{.line_color = color_, .line_width = line_width_}));
      layer->draw(lines_.back());
    }
  }
}

} // namespace cairowindow::draw
