#include "sys.h"
#include "LayerRegion.h"
#include "Layer.h"
#include "debug.h"

namespace cairowindow {

void LayerRegion::draw(Layer* layer)
{
  DoutEntering(dc::notice, "LayerRegion::draw(" << layer << ") [" << this << "]");

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
  layer->add_area(stroke_extents_.area());
  layer->window_update(stroke_extents_);
}

LayerRegion::~LayerRegion()
{
  DoutEntering(dc::notice, "~LayerRegion() [" << this << "]");
  layer_->remove(this);
  layer_->window_update(stroke_extents_);
}

} // namespace cairowindow
