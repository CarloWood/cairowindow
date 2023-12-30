#include "sys.h"
#include "PlotArea.h"
#include "cairowindow/Range.h"
#include "cairowindow/Layer.h"
#include "utils/almost_equal.h"
#include <algorithm>
#include "debug.h"

namespace cairowindow::draw {

// Draw a line from x1,y1 to x2,y2 and draw k ticks perpendicular to that line, pointing clockwise.
void PlotArea::draw_axis(cairo_t* cr, double x1, double y1, double x2, double y2, int k)
{
  cairo_move_to(cr, x1, y1);
  cairo_line_to(cr, x2, y2);
  double sx = (x2 - x1) / k;
  double sy = (y2 - y1) / k;
  double sl = std::max(std::abs(sx), std::abs(sy));
  double tx = -sy * tick_length_ / sl;
  double ty = sx * tick_length_ / sl;
  // Draw the ticks.
  for (int t = 1; t < k; ++t)
  {
    x1 += sx;
    y1 += sy;
    cairo_move_to(cr, x1, y1);
    cairo_line_to(cr, x1 + tx, y1 + ty);
  }
}

//static
int PlotArea::calculate_range_ticks(plot::Range& range)
{
  double diff = range.max() - range.min();
  double order = std::floor(std::log10(diff) + 1e-6);
  double spacing = std::pow(10.0, order);
  double ticks = diff / spacing;
  if (ticks < 1.00001)
    spacing *= 0.1;
  else if (ticks < 2.00001)
    spacing *= 0.2;
  else if (ticks < 5.00001)
    spacing *= 0.5;
  range.round_to(spacing);
  return std::round((range.max() - range.min()) / spacing);
}

void PlotArea::draw_regions_on(Layer* layer)
{
  std::array<int, number_of_axes> k;
  for (int axis = 0; axis < 2; ++axis)
  {
    plot::Range range{range_[axis][min_range], range_[axis][max_range]};
    k[axis] = calculate_range_ticks(range);
  }

  if (draw_grid_)
  {
    grid_.set_ticks(k);
    layer->draw(&grid_);
  }

  auto x_axis_min = [this, k = k[x_axis]](cairo_t* cr) -> StrokeExtents
  {
    set_line_style(cr);
    draw_axis(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y() + geometry_.height(),
        geometry_.offset_x(), geometry_.offset_y() + geometry_.height(), k);
    return stroke(cr);
  };

  auto x_axis_max = [this, k = k[x_axis]](cairo_t* cr) -> StrokeExtents
  {
    set_line_style(cr);
    draw_axis(cr, geometry_.offset_x(), geometry_.offset_y(),
        geometry_.offset_x() + geometry_.width(), geometry_.offset_y(), k);
    return stroke(cr);
  };

  auto y_axis_min = [this, k = k[y_axis]](cairo_t* cr) -> StrokeExtents
  {
    set_line_style(cr);
    draw_axis(cr, geometry_.offset_x(), geometry_.offset_y() + geometry_.height(),
        geometry_.offset_x(), geometry_.offset_y(), k);
    return stroke(cr);
  };

  auto y_axis_max = [this, k = k[y_axis]](cairo_t* cr) -> StrokeExtents
  {
    set_line_style(cr);
    draw_axis(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y(),
        geometry_.offset_x() + geometry_.width(), geometry_.offset_y() + geometry_.height(), k);
    return stroke(cr);
  };

  axes_[x_axis][min_range] = std::make_shared<LayerRegion>(x_axis_min);
  axes_[x_axis][max_range] = std::make_shared<LayerRegion>(x_axis_max);
  axes_[y_axis][min_range] = std::make_shared<LayerRegion>(y_axis_min);
  axes_[y_axis][max_range] = std::make_shared<LayerRegion>(y_axis_max);

  layer->draw(axes_[x_axis][min_range]);
  layer->draw(axes_[x_axis][max_range]);
  layer->draw(axes_[y_axis][min_range]);
  layer->draw(axes_[y_axis][max_range]);
}

void PlotArea::set_range(int axis, double range_min, double range_max)
{
  range_[axis][min_range] = range_min;
  range_[axis][max_range] = range_max;
}

} // namespace cairowindow::draw
