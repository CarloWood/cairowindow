#pragma once

#include "Line.h"
#include "MultiRegion.h"
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
  cairowindow::Rectangle geometry_;     // The geometry passed to the constructor.
                                        // This is the path used for the large rectangle around the plot area.
  std::array<int, number_of_axes> ticks_;
  std::vector<std::shared_ptr<Line>> lines_;

  // Implementation of MultiRegion.
  void draw_regions_on(Layer* layer) override;

 public:
  Grid(cairowindow::Rectangle const& geometry, GridStyle style) :
    MultiRegion(style.color, style.line_width), geometry_(geometry) { }

  void set_ticks(int axis, int ticks)
  {
    ticks_[axis] = ticks;
  }

  void set_geometry(cairowindow::Rectangle const& geometry)
  {
    geometry_ = geometry;
  }

  std::array<int, number_of_axes> const& ticks() const
  {
    return ticks_;
  }
};

} // namespace cairowindow::draw
