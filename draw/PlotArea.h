#pragma once

#include "cairowindow/MultiRegion.h"
#include "cairowindow/LayerRegion.h"
#include "cairowindow/Rectangle.h"
#include "cairowindow/Color.h"
#include <array>

namespace cairowindow::plot {

constexpr int x_axis = 0;
constexpr int y_axis = 1;

constexpr int min_range = 0;
constexpr int max_range = 1;

class Range;

} // namespace cairowindow::plot

namespace cairowindow::draw {

struct PlotAreaStyle
{
  Color axes_color = color::black;
  double axes_line_width = 1.0;
};

class PlotArea : public MultiRegion
{
  // Lets not get confused with draw::Rectangle (in case that is #include-d).
  using Rectangle = cairowindow::Rectangle;

 public:
  static constexpr int x_axis = plot::x_axis;
  static constexpr int y_axis = plot::y_axis;
  static constexpr int min_range = plot::min_range;
  static constexpr int max_range = plot::max_range;
  static constexpr int number_of_axes = 2;

 private:
  Rectangle geometry_;          // The geometry passed to the constructor. This is the path used for the large rectangle around the plot area.
  double tick_length_;

  std::array<std::array<LayerRegion, 2>, number_of_axes> axes_;
  std::array<std::array<double, 2>, number_of_axes> range_{{{0, 1}, {0, 1}}};

 public:
  PlotArea(Rectangle const& geometry, PlotAreaStyle style) :
    MultiRegion(style.axes_color, style.axes_line_width), geometry_(geometry), tick_length_(geometry.width() / 100.0) { }

  void set_range(int axis, double range_min, double range_max);

  Rectangle const& geometry() const { return geometry_; }

  static int calculate_range_ticks(plot::Range& range);

 private:
  void draw_axis(cairo_t* cr, double x1, double y1, double x2, double y2, int k);
  void draw_regions_on(Layer* layer) override;
};

} // namespace cairowindow::draw
