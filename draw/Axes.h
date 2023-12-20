#pragma once

#include "cairowindow/MultiRegion.h"
#include "cairowindow/LayerRegion.h"
#include "cairowindow/Rectangle.h"
#include "cairowindow/Color.h"

namespace cairowindow::draw {

class Axes : public MultiRegion
{
  // Lets not get confused with draw::Rectangle (in case that is #include-d).
  using Rectangle = cairowindow::Rectangle;

 private:
  Rectangle geometry_;          // The geometry passed to the constructor. This is the path that is used for the large rectangle of the axes.
  LayerRegion region_;

 public:
  Axes(Rectangle const& geometry, Color const& color, double line_width = 1.0) :
    MultiRegion(color, line_width), geometry_(geometry) { }

 private:
  void draw_regions_on(Layer* layer) override;
};

} // namespace cairowindow::draw
