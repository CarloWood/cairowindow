#include "sys.h"
#include "LayerRegion.h"
#include "Layer.h"
#include "debug.h"
#ifdef CWDEBUG
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow {

void LayerRegion::draw(Layer* layer)
{
  DoutEntering(dc::cairowindow, "LayerRegion::draw(" << layer << ") [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif

  // Only call draw once.
  ASSERT(!layer_);
  layer_ = layer;

  cairo_t* cr = layer->cr();
  cairo_save(cr);
  // Apply layer offset, if any.
  cairo_translate(cr, -layer->offset_x(), -layer->offset_y());

  // For now only a single user_draw() call per region is allowed.
  // The reason for this is that a region is supposed to a single rectangle, so doing
  // multiple draw calls could accidently cause large rectangles with a lot of empty space.
  ASSERT(!stroke_extents_.is_defined());
  // Call the drawing function of the user.
  stroke_extents_ = std::move(redraw(cr));

  cairo_restore(cr);

  if (stroke_extents_.is_defined())
  {
    layer->add_area(stroke_extents_.area());
    layer->window_update(stroke_extents_);
  }
}

StrokeExtents LayerRegion::redraw(cairo_t* cr)
{
  DoutEntering(dc::cairowindow, "LayerRegion::redraw(cr) [" << this << "]");

  StrokeExtents result = draw_ ? draw_(cr) : do_draw(cr);
  if (!result.clip(layer_->geometry()))
    return {};
  return result;
}

LayerRegion::~LayerRegion()
{
  DoutEntering(dc::cairowindow, "~LayerRegion() [" << this << "]");
  if (layer_)
  {
    layer_->remove(this);
    layer_->window_update(stroke_extents_);
  }
}

} // namespace cairowindow
