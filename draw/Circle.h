#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"

namespace cairowindow::draw {

struct CircleStyle
{
  // These must correspond 1:1 with the first four elements of ShapeStyle!
  Color line_color = color::black;
  Color fill_color = color::transparent;
  double line_width = 1.0;
  bool at_corner = true;                        // The center of the circle is normally at the top-left corner.
};

class Circle : public Shape
{
 public:
  Circle(Rectangle const& geometry, CircleStyle style) :
    Shape(geometry, ShapeStyle{style.line_color, style.fill_color, style.line_width, style.at_corner, ellipse}) { }

  CircleStyle& style() { return reinterpret_cast<CircleStyle&>(style_); }
  CircleStyle const& style() const { return reinterpret_cast<CircleStyle const&>(style_); }
};

} // namespace cairowindow::draw
