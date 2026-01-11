#pragma once

#include "cairowindow/Rectangle.h"
#include <memory>

namespace cairowindow::draw {
// Forward declare the draw object.
class Rectangle;
} // namespace cairowindow::draw

namespace cairowindow::plot {
// Forward declaration.
class Plot;

//-----------------------------------------------------------------------------
// plot::Rectangle
//
// A handle keeping a plotted Rectangle alive.
// Returned by Plot::create_rectangle(layer, rectangle_style, <args to construct a plot::Rectangle>).
//
class Rectangle : public cairowindow::Rectangle
{
 public:
  explicit Rectangle(cairowindow::Rectangle const& rectangle) : cairowindow::Rectangle(rectangle) { }
  using cairowindow::Rectangle::Rectangle;

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Rectangle> draw_object_;
};

//
//-----------------------------------------------------------------------------

} // namespace cairowindow::plot
