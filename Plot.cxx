#include "sys.h"
#include "Layer.h"
#include "Plot.h"

namespace cairowindow::plot {

// Calculate the axes geometry from the full plot geometry.
Rectangle Plot::axes_geometry(Rectangle const& geometry)
{
  // Use fixed-size margins for now.
  constexpr int top_margin = 50;
  constexpr int bottom_margin = 60;
  constexpr int left_margin = 100;
  constexpr int right_margin = 20;

  return { geometry.offset_x() + left_margin, geometry.offset_y() + top_margin,
           geometry.width() - left_margin - right_margin, geometry.height() - top_margin - bottom_margin };
}

void Plot::add_to(boost::intrusive_ptr<Layer> const& layer)
{
  layer->draw(&plot_area_);
  if (title_.is_defined())
    layer->draw(&title_);
  if (xlabel_.is_defined())
    layer->draw(&xlabel_);
  if (ylabel_.is_defined())
    layer->draw(&ylabel_);
}

} // namespace cairowindow::plot
