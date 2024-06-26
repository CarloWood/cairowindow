#pragma once

#include "Grid.h"
#include "cairowindow/LayerRegion.h"
#include "cairowindow/Rectangle.h"
#include "cairowindow/Color.h"
#include <array>

namespace cairowindow {
class Range;
namespace plot {

constexpr int x_axis = 0;
constexpr int y_axis = 1;

constexpr int min_range = 0;
constexpr int max_range = 1;

} // namespace plot
} // namespace cairowindow

namespace cairowindow::draw {

struct PlotAreaStyle
{
  Color axes_color = color::black;
  double axes_line_width = 1.0;
  GridStyle grid{};
};

class PlotArea : public MultiRegion
{
 public:
  static constexpr int x_axis = plot::x_axis;
  static constexpr int y_axis = plot::y_axis;
  static constexpr int min_range = plot::min_range;
  static constexpr int max_range = plot::max_range;
  static constexpr int number_of_axes = Grid::number_of_axes;

 private:
  cairowindow::Rectangle geometry_;             // The geometry passed to the constructor.
                                                // This is the path used for the large rectangle around the plot area.
  double tick_length_;
  bool draw_grid_;
  draw::Grid grid_;

  std::array<std::array<std::shared_ptr<LayerRegion>, 2>, number_of_axes> axes_;
  std::array<std::array<double, 2>, number_of_axes> range_{{{0, 1}, {0, 1}}};

 public:
  PlotArea(cairowindow::Rectangle const& geometry, PlotAreaStyle style) :
    MultiRegion(style.axes_color, style.axes_line_width), geometry_(geometry), tick_length_(geometry.width() / 100.0),
    draw_grid_(!style.grid.color.is_transparent()), grid_(geometry, style.grid) { }

  void set_range(int axis, double range_min, double range_max, int ticks);
  void set_geometry(cairowindow::Rectangle const& geometry)
  {
    geometry_ = geometry;
    grid_.set_geometry(geometry);
  }

  cairowindow::Rectangle const& geometry() const { return geometry_; }

  static int calculate_range_ticks(Range& range);

 private:
  void draw_axis(cairo_t* cr, double x1, double y1, double x2, double y2, int k);
  void draw_regions_on(Layer* layer) override;

 public:
#ifdef CWDEBUG
  friend std::ostream& operator<<(std::ostream& os, PlotArea const* plot_area_ptr)
  {
    os << "PlotArea*";
    return os;
  }
#endif
};

} // namespace cairowindow::draw
