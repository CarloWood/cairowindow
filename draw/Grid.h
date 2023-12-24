#pragma once

#include "Line.h"
#include "cairowindow/MultiRegion.h"
#include <array>
#include <vector>

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
  std::vector<std::unique_ptr<Line>> lines_;

  // Implementation of MultiRegion.
  void draw_regions_on(Layer* layer) override;

 public:
  Grid(Rectangle const& geometry, GridStyle style) :
    MultiRegion(style.color, style.line_width), geometry_(geometry) { }

  void set_ticks(std::array<int, number_of_axes> const& k)
  {
    ticks_ = k;
  }
};

} // namespace cairowindow::draw
