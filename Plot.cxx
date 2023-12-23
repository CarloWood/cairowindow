#include "sys.h"
#include "Layer.h"
#include "Plot.h"
#include <cmath>

namespace cairowindow::plot {

// Calculate the axes geometry from the full plot geometry.
Rectangle Plot::axes_geometry(Rectangle const& geometry, double axes_line_width)
{
  // Use fixed-size margins for now.
  constexpr int top_margin = 50;
  constexpr int bottom_margin = 60;
  constexpr int left_margin = 100;
  constexpr int right_margin = 20;

  double top_left_x = geometry.offset_x() + left_margin;
  double top_left_y = geometry.offset_y() + top_margin;
  double bottom_right_x = geometry.offset_x() + geometry.width() - right_margin;
  double bottom_right_y = geometry.offset_y() + geometry.height() - bottom_margin;

  int lwi = std::round(axes_line_width);
  if (lwi == axes_line_width)
  {
    if ((lwi & 1) == 0)
    {
      top_left_x = std::round(top_left_x);
      top_left_y = std::round(top_left_y);
      bottom_right_x = std::round(bottom_right_x);
      bottom_right_y = std::round(bottom_right_y);
    }
    else
    {
      top_left_x = std::ceil(top_left_x) - 0.5;
      top_left_y = std::ceil(top_left_y) - 0.5;
      bottom_right_x = std::floor(bottom_right_x) + 0.5;
      bottom_right_y = std::floor(bottom_right_y) + 0.5;
    }
  }

  return { top_left_x, top_left_y, bottom_right_x - top_left_x, bottom_right_y - top_left_y };
}

void Plot::add_to(boost::intrusive_ptr<Layer> const& layer)
{
  // Set a title.
  if (title_.is_defined())
    layer->draw(&title_);

  // Set ranges on the plot area and draw it.
  for (int axis = x_axis; axis <= y_axis; ++axis)
    plot_area_.set_range(axis, range_[axis].min(), range_[axis].max());
  layer->draw(&plot_area_);

  // Draw axis labels.
  double max_height = 0;
  double max_width = 0;
  for (int axis = x_axis; axis <= y_axis; ++axis)
  {
    double label = range_[axis].min();
    double delta = range_[axis].size() / range_ticks_[axis];
    int precision = -std::floor(std::log10(delta) + 1e-6);
    double x = plot_area_.geometry().offset_x();
    double y = plot_area_.geometry().offset_y() + plot_area_.geometry().height();
    if (axis == x_axis)
      y -= XLabelStyleDefaults::offset;
    else
      x += XLabelStyleDefaults::offset;
    for (int tick = 0; tick <= range_ticks_[axis]; ++tick)
    {
      std::ostringstream label_str;
      if (precision > 0)
        label_str << std::fixed << std::setprecision(precision);
      label_str << label;
      labels[axis].emplace_back(std::make_unique<draw::Text>(label_str.str(),
            x, y, LabelStyle{.position = (axis == x_axis) ? draw::centered_below : draw::centered_left_of}));
      layer->draw(labels[axis].back());
      StrokeExtents label_extents = labels[axis].back()->stroke_extents();
      if (axis == x_axis)
        max_height = std::max(max_height, label_extents.height());
      else
        max_width = std::max(max_width, label_extents.width());

      label += delta;
      if (axis == x_axis)
        x += plot_area_.geometry().width() / range_ticks_[axis];
      else
        y -= plot_area_.geometry().height() / range_ticks_[axis];
    }
  }

  if (xlabel_.is_defined())
  {
    xlabel_.move_pos(0, max_height);
    layer->draw(&xlabel_);
  }

  if (ylabel_.is_defined())
  {
    ylabel_.move_pos(-max_width, 0);
    layer->draw(&ylabel_);
  }
}

} // namespace cairowindow::plot
