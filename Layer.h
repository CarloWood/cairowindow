#pragma once

#include "Color.h"
#include "Rectangle.h"
#include <cairo/cairo.h>
#include <functional>

namespace cairowindow {

class Window;

class Layer
{
 private:
  Window* window_;
  cairo_surface_t* drawing_surface_;
  cairo_t* drawing_cr_;
  Rectangle rectangle_;

 public:
  Layer(cairo_surface_t* x11_surface, Rectangle const& rectangle, cairo_content_t content, Color color, Window* window);
  ~Layer();

  void draw(std::function<Rectangle(cairo_t*)> user_draw);

  cairo_surface_t* drawing_surface() const { return drawing_surface_; }
  double offset_x() const { return rectangle_.offset_x(); }
  double offset_y() const { return rectangle_.offset_y(); }
};

} // namespace cairowindow
