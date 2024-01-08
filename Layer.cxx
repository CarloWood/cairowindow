#include "sys.h"
#include "Layer.h"
#include "Window.h"
#include "debug.h"
#ifdef CWDEBUG
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow {

Layer::Layer(cairo_surface_t* x11_surface, Rectangle const& rectangle, cairo_content_t content, Color color, Window* window
    COMMA_DEBUG_ONLY(std::string debug_name)) :
  color_(color), window_(window), geometry_(rectangle), region_areas_(0.0)
{
  DoutEntering(dc::cairowindow, "Layer::Layer(" << x11_surface << ", " << rectangle << ", " << content << ", " << color << ", " << window <<
      ", \"" << debug_name << "\") [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  // Create an off-screen surface for tripple buffering.
  surface_ = cairo_surface_create_similar(x11_surface, content, rectangle.width(), rectangle.height()
      COMMA_CWDEBUG_ONLY("Layer::surface_:\"" + debug_name + "\""));
  cr_ = cairo_create(surface_
      COMMA_CWDEBUG_ONLY("Layer::cr_:\"" + debug_name + "\""));
#ifdef CAIROWINDOW_DEBUGWINDOW
  debug_window_.start(surface_, rectangle.width(), rectangle.height(), debug_name);
#endif

  if (content == CAIRO_CONTENT_COLOR)
    color_.set_opaque();

  // Set background color for off-screen surface.
  cairo_set_source_rgba(cr_, color_.red(), color_.green(), color_.blue(), color_.alpha());

  // Fill surface with background color.
  cairo_paint(cr_);
}

Layer::~Layer()
{
  DoutEntering(dc::cairowindow, "Layer::~Layer() [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
#ifdef CAIROWINDOW_DEBUGWINDOW
  debug_window_.terminate();
#endif
  // Clean up.
  cairo_destroy(cr_);
  cairo_surface_destroy(surface_);
}

void Layer::redraw(cairo_t* cr, StrokeExtents const& stroke_extents)
{
  DoutEntering(dc::cairowindow, "Layer::redraw(" << cr << ", " << stroke_extents << ") [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  bool clip_rectangle_set = false;
  auto last = regions_.end();
  for (auto it = regions_.begin(); it != last;)
  {
    auto layer_region = it->lock();
    if (!layer_region)
    {
      if (it == --last)
        break;
      std::iter_swap(it, last);
    }
    else
    {
      if (stroke_extents.overlaps(layer_region->stroke_extents()))
      {
        if (!clip_rectangle_set)
        {
          cairo_save(cr);
          stroke_extents.set_path(cr);
          cairo_clip(cr);
          clip_rectangle_set = true;
        }
        cairo_save(cr);
        layer_region->redraw(cr);
        cairo_restore(cr);
      }
      ++it;
    }
  }
  if (clip_rectangle_set)
    cairo_restore(cr);
  regions_.erase(last, regions_.end());
}

void Layer::remove(LayerRegion const* layer_region)
{
  DoutEntering(dc::cairowindow, "Layer::remove(" << layer_region << ") [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  regions_.erase(std::remove_if(regions_.begin(), regions_.end(),
        [](std::weak_ptr<LayerRegion> const& wp){ return wp.expired(); }), regions_.end());
  // Replace rectangle with background color.
  cairo_save(cr_);
  cairo_translate(cr_, -geometry_.offset_x(), -geometry_.offset_y());
  cairo_set_operator(cr_, CAIRO_OPERATOR_SOURCE);
  cairo_set_source_rgba(cr_, color_.red(), color_.green(), color_.blue(), color_.alpha());
  StrokeExtents const& stroke_extents = layer_region->stroke_extents();
  stroke_extents.set_path(cr_);
  cairo_fill(cr_);
  cairo_set_operator(cr_, CAIRO_OPERATOR_OVER);
  redraw(cr_, stroke_extents);
  cairo_restore(cr_);
}

} // namespace cairowindow
