#pragma once

#include "Color.h"
#include "Rectangle.h"
#include "LayerRegion.h"
#include "Window.h"
#include "utils/AIRefCount.h"
#include <cairo/cairo.h>
#include <vector>
#include "debug.h"

namespace cairowindow {

template<typename Type>
concept LayerRegionType = std::is_base_of_v<LayerRegion, Type>;

class Layer : public AIRefCount
{
 private:
  Color color_;
  Window* window_;
  cairo_surface_t* surface_;
  cairo_t* cr_;
  Rectangle rectangle_;
  std::vector<LayerRegion*> regions_;
  double region_areas_; // Total area of all regions_.
#ifdef CAIROWINDOW_DEBUGWINDOW
  DebugWindow debug_window_;
#endif

 public:
  Layer(cairo_surface_t* x11_surface, Rectangle const& rectangle, cairo_content_t content, Color color, Window* window
      COMMA_DEBUG_ONLY(std::string debug_name));
  ~Layer();

  Layer(Layer const& layer) = delete;

  template<LayerRegionType LRT, typename... ARGS>
  boost::intrusive_ptr<LRT> create_layer_region(ARGS&&... args)
  {
    DoutEntering(dc::notice, "Layer::create_layer_region<" <<
        libcwd::type_info_of<LRT>().demangled_name() << ">(" << join(", ", args...) << ") [" << this << "]");
    boost::intrusive_ptr<LRT> region = new LRT(this, std::forward<ARGS>(args)...);
    regions_.push_back(region.get());
    return region;
  }

  void remove(LayerRegion* region);

  void window_update(Rectangle const& rectangle)
  {
#ifdef CAIROWINDOW_DEBUGWINDOW
    debug_window_.update(rectangle);
#endif
    window_->update(rectangle);
  }

  void redraw(cairo_t* cr, Rectangle const& rectangle);
  void add_area(double area) { region_areas_ += area; }

  cairo_surface_t* surface() const { return surface_; }
  cairo_t* cr() const { return cr_; }
  double offset_x() const { return rectangle_.offset_x(); }
  double offset_y() const { return rectangle_.offset_y(); }
  double area() const { return region_areas_; }
};

} // namespace cairowindow
