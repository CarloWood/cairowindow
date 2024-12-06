#include "sys.h"
#include "utils/generate_unique_filename.h"
#include "Printable.h"
#include "Layer.h"
#include <cairo/cairo-svg.h>
#ifdef CWDEBUG
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow {

Printable::~Printable()
{
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  if (svg_surface_)
  {
    cairo_surface_finish(svg_surface_);
    cairo_destroy(svg_cr_);
    cairo_surface_destroy(svg_surface_);
  }
}

void Printable::create_svg_surface(std::string svg_filename, bool overwrite COMMA_CWDEBUG_ONLY(std::string debug_name))
{
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  if (!overwrite)
  {
    // Note, if you get
    // error: no member named 'generate_unique_filename' in namespace 'utils'
    // Add 'find_package(Boost COMPONENTS filesystem)' to the CMakeLists.txt file of the root project.
    svg_filename = utils::generate_unique_filename(svg_filename);
  }
  Rectangle const& geometry = this->geometry();
  svg_surface_ = cairo_svg_surface_create(svg_filename.c_str(), geometry.width(), geometry.height() COMMA_CWDEBUG_ONLY(debug_name));
  svg_cr_ = cairo_create(svg_surface_ COMMA_CWDEBUG_ONLY("Plot::svg_cr_:\"" + debug_name + "\""));
  // Define a clip path for the whole area.
  cairo_rectangle(svg_cr_, 0, 0, geometry.width(), geometry.height());
  cairo_clip(svg_cr_);
  // And fill it entirely with a white background.
  cairo_set_source_rgb(svg_cr_, 1.0, 1.0, 1.0);
  cairo_paint(svg_cr_);
  // Set translate to move all drawing to the top-left corner of this SVG surface.
  cairo_translate(svg_cr_, -geometry.offset_x(), -geometry.offset_y());
}

void Printable::draw_layer_region_on(boost::intrusive_ptr<Layer> const& layer, std::shared_ptr<LayerRegion> const& layer_region)
{
  if (need_print_)
  {
    ASSERT(svg_cr_);
    layer->start_printing_to(svg_cr_);
  }
  layer->draw(layer_region);
  if (need_print_)
    layer->stop_printing();
}

void Printable::draw_multi_region_on(boost::intrusive_ptr<Layer> const& layer, draw::MultiRegion* multi_region)
{
  if (need_print_)
  {
    ASSERT(svg_cr_);
    layer->start_printing_to(svg_cr_);
  }
  layer->draw(multi_region);
  if (need_print_)
    layer->stop_printing();
}

} // namespace cairowindow
