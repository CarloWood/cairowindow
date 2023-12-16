#include "sys.h"
#include "Layer.h"
#include "Window.h"
#include "debug.h"

namespace cairowindow {

Layer::Layer(cairo_surface_t* x11_surface, Rectangle const& rectangle, cairo_content_t content, Color color, Window* window) :
  window_(window), rectangle_(rectangle), region_areas_(0.0)
{
  // Create an off-screen surface for tripple buffering.
  surface_ = cairo_surface_create_similar(x11_surface, content, rectangle.width(), rectangle.height());
  cr_ = cairo_create(surface_);

  if (content == CAIRO_CONTENT_COLOR)
  {
    // Set background color for off-screen surface.
    cairo_set_source_rgb(cr_, color.red(), color.green(), color.blue());

    // Fill surface with background color.
    cairo_paint(cr_);
  }
  else if (color.alpha() > 0)
  {
    // Set background color for off-screen surface.
    cairo_set_source_rgba(cr_, color.red(), color.green(), color.blue(), color.alpha());

    // Fill surface with background color.
    cairo_paint(cr_);
  }
}

Layer::~Layer()
{
  // Clean up.
  cairo_destroy(cr_);
  cairo_surface_destroy(surface_);
}

void Layer::redraw(cairo_t* cr, Rectangle const& rectangle)
{
  for (auto const& region_ptr : regions_)
    if (rectangle.overlaps(region_ptr->rectangle()))
      region_ptr->redraw(cr);
}

} // namespace cairowindow
