#pragma once

#include "Line.h"
#include "MultiRegion.h"
#include <array>
#include <vector>
#include <memory>

namespace cairowindow::draw {

#define cairowindow_Grid_FOREACH_MEMBER(X, ...) \
  X(Color, color, Color{}, __VA_ARGS__) \
  X(double, line_width, -1.0, __VA_ARGS__)

#define cairowindow_Grid_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Grid_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for GridStyle.
struct GridStyleParamsDefault
{
  static constexpr Color color = color::transparent;
  static constexpr double line_width = 1.0;
};

// Declare GridStyle.
DECLARE_STYLE(Grid, GridStyleParamsDefault);

class Grid : public MultiRegion
{
 public:
 static constexpr int number_of_axes = 2;

 private:
  cairowindow::Geometry geometry_;     // The geometry passed to the constructor.
                                        // This is the path used for the large rectangle around the plot area.
  std::array<int, number_of_axes> ticks_;
  std::vector<std::shared_ptr<Line>> lines_;

  // Implementation of MultiRegion.
  void draw_regions_on(Layer* layer) override;

 public:
  Grid(cairowindow::Geometry const& geometry, GridStyle style) :
    MultiRegion(style.color(), style.line_width()), geometry_(geometry) { }

  void set_ticks(int axis, int ticks)
  {
    ticks_[axis] = ticks;
  }

  void set_geometry(cairowindow::Geometry const& geometry)
  {
    geometry_ = geometry;
  }

  std::array<int, number_of_axes> const& ticks() const
  {
    return ticks_;
  }
};

} // namespace cairowindow::draw
