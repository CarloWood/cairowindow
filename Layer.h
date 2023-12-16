#pragma once

#include "Color.h"
#include "Rectangle.h"
#include "LayerRegion.h"
#include "Window.h"
#include "utils/AIRefCount.h"
#include <cairo/cairo.h>
#include <vector>

namespace cairowindow {

template<typename Type>
concept LayerRegionType = std::is_base_of_v<LayerRegion, Type>;

class Layer : public AIRefCount
{
 private:
  Window* window_;
  cairo_surface_t* surface_;
  cairo_t* cr_;
  Rectangle rectangle_;
  std::vector<boost::intrusive_ptr<LayerRegion>> regions_;
  double region_areas_; // Total area of all regions_.

 public:
  Layer(cairo_surface_t* x11_surface, Rectangle const& rectangle, cairo_content_t content, Color color, Window* window);
  ~Layer();

  Layer(Layer const& layer) = delete;

  template<LayerRegionType LRT, typename... ARGS>
  boost::intrusive_ptr<LRT> create_layer_region(ARGS&&... args)
  {
    boost::intrusive_ptr<LRT> region = new LRT(this, std::forward<ARGS>(args)...);
    regions_.push_back(region);
    return region;
  }

  void window_update(Rectangle const& rectangle) const { window_->update(rectangle); }

  void redraw(cairo_t* cr, Rectangle const& rectangle);
  void add_area(double area) { region_areas_ += area; }

  cairo_surface_t* surface() const { return surface_; }
  cairo_t* cr() const { return cr_; }
  double offset_x() const { return rectangle_.offset_x(); }
  double offset_y() const { return rectangle_.offset_y(); }
  double area() const { return region_areas_; }
};

} // namespace cairowindow
