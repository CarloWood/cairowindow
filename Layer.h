#pragma once

#include "Color.h"
#include "Rectangle.h"
#include "LayerRegion.h"
#include "MultiRegion.h"
#include "Window.h"
#include "utils/AIRefCount.h"
#include <cairo/cairo.h>
#include <vector>
#include <memory>
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
  Rectangle geometry_;
  std::vector<std::weak_ptr<LayerRegion>> regions_;
  double region_areas_; // Total area of all regions_.
#ifdef CAIROWINDOW_DEBUGWINDOW
  DebugWindow debug_window_;
#endif

 public:
  Layer(cairo_surface_t* x11_surface, Rectangle const& rectangle, cairo_content_t content, Color color, Window* window
      COMMA_DEBUG_ONLY(std::string debug_name));
  ~Layer();

  Layer(Layer const& layer) = delete;

  void draw(std::shared_ptr<LayerRegion> const& layer_region)
  {
    regions_.emplace_back(layer_region);
    layer_region->draw(this);
  }

  void draw(MultiRegion* multi_region)
  {
    multi_region->draw_regions_on(this);
  }

  void remove(LayerRegion const* layer_region);

  void window_update(StrokeExtents const& stroke_extents)
  {
#ifdef CAIROWINDOW_DEBUGWINDOW
    debug_window_.update(stroke_extents);
#endif
    window_->update(stroke_extents);
  }

  void redraw(cairo_t* cr, StrokeExtents const& stroke_extents);
  void add_area(double area) { region_areas_ += area; }

  cairo_surface_t* surface() const { return surface_; }
  cairo_t* cr() const { return cr_; }
  double offset_x() const { return geometry_.offset_x(); }
  double offset_y() const { return geometry_.offset_y(); }
  double area() const { return region_areas_; }
  // Return the geometry of the layer.
  Rectangle const& geometry() const { return geometry_; }

#ifdef CWDEBUG
  friend std::ostream& operator<<(std::ostream& os, Layer const* layer_ptr)
  {
    os << "Layer*";
    return os;
  }
#endif
};

} // namespace cairowindow
