#include "sys.h"
#include "PlotArea.h"
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

void PlotArea::draw_regions_on(Layer* layer)
{
  std::array<int, number_of_axes> k;
  for (int axis = 0; axis < 2; ++axis)
  {
    double range = range_[axis][max_range] - range_[axis][min_range];
    double order = std::floor(std::log10(range));
    double spacing = std::pow(10.0, order);
    double ticks = range / spacing;
    if (ticks < 2.0)
      spacing *= 0.2;
    else if (ticks < 5.0)
      spacing *= 0.5;
    // The range that is set must already be preprocessed to be an integer times the spacing that we calculated here.
    ASSERT(utils::almost_equal(std::round(range_[axis][min_range] / spacing), range_[axis][min_range] / spacing, 10e-4));
    ASSERT(utils::almost_equal(std::round(range_[axis][max_range] / spacing), range_[axis][max_range] / spacing, 10e-4));
    k[axis] = std::round(range / spacing);
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

  axes_[x_axis][min_range].draw(layer, x_axis_min);
  axes_[x_axis][max_range].draw(layer, x_axis_max);
  axes_[y_axis][min_range].draw(layer, y_axis_min);
  axes_[y_axis][max_range].draw(layer, y_axis_max);
}

void PlotArea::set_range(int axis, double range_min, double range_max)
{
  range_[axis][min_range] = range_min;
  range_[axis][max_range] = range_max;
}

} // namespace cairowindow::draw
