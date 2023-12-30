#pragma once

#include "Line.h"
#include "cairowindow/MultiRegion.h"
#include <array>
#include <vector>
#include <memory>

namespace cairowindow::draw {

struct GridStyle
{
  Color color = color::transparent;
  double line_width = 1.0;
};

class Grid : public MultiRegion
{
 public:
  static constexpr int number_of_axes = 2;

 private:
  Rectangle geometry_;          // The geometry passed to the constructor. This is the path used for the large rectangle around the plot area.
  std::array<int, number_of_axes> ticks_;
  std::vector<std::shared_ptr<Line>> lines_;

  // Implementation of MultiRegion.
  void draw_regions_on(Layer* layer) override;

 public:
  Grid(Rectangle const& geometry, GridStyle style) :
    MultiRegion(style.color, style.line_width), geometry_(geometry) { }

  void set_ticks(std::array<int, number_of_axes> const& k)
  {
    ticks_ = k;
  }

  void set_geometry(Rectangle const& geometry)
  {
    geometry_ = geometry;
  }
};

} // namespace cairowindow::draw
