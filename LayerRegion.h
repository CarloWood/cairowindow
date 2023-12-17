#pragma once

#include "Rectangle.h"
#include "utils/AIRefCount.h"
#include <cairo/cairo.h>
#include <functional>

namespace cairowindow {

class Layer;

class LayerRegion : public AIRefCount
{
 private:
  Layer* layer_;                                // The Layer that this is a region of. Drawing happens to the surface of this layer.
  Rectangle rectangle_;                         // The region area (as returned by draw).
  std::function<Rectangle(cairo_t*)> draw_;     // A copy the argument of draw; used to redraw the region when necessary.

 private:
  virtual Rectangle do_draw(cairo_t* cr) { Dout(dc::warning, "Calling unimplemented do_draw()"); return {}; }

 public:
  LayerRegion(Layer* layer) : layer_(layer) { }

  void draw();
  void draw(std::function<Rectangle(cairo_t*)> user_draw)
  {
    draw_ = user_draw;
    draw();
  }

  Rectangle redraw(cairo_t* cr)
  {
    if (draw_)
      return draw_(cr);

    return do_draw(cr);
  }

  Rectangle const& rectangle() const { return rectangle_; }
};

} // namespace cairowindow
