#pragma once

#include <cairo/cairo.h>
#include <boost/intrusive_ptr.hpp>
#include <memory>
#include <string>
#include "debug.h"

namespace cairowindow {

class Layer;
class LayerRegion;
class Rectangle;
namespace draw {
class MultiRegion;
} // namespace draw

class Printable
{
 protected:
  cairo_surface_t* svg_surface_ = nullptr;
  cairo_t* svg_cr_ = nullptr;
  bool need_print_{false};

 public:
  virtual ~Printable();

  void create_svg_surface(std::string svg_filename, bool overwrite COMMA_CWDEBUG_ONLY(std::string debug_name));

  cairo_t* svg_cr() const { return svg_cr_; }

  void set_need_print() { need_print_ = true; }
  void reset_need_print() { need_print_ = false; }
  bool need_print() const { return need_print_; }

  // Calls layer->draw(layer_region), but also handles printing.
  void draw_layer_region_on(boost::intrusive_ptr<Layer> const& layer, std::shared_ptr<LayerRegion> const& layer_region);
  void draw_multi_region_on(boost::intrusive_ptr<Layer> const& layer, draw::MultiRegion* multi_region);

  virtual Rectangle const& geometry() const = 0;
};

} // namespace cairowindow
