#include "sys.h"
#include "Layer.h"
#include "Window.h"

namespace cairowindow {

Layer::Layer(cairo_surface_t* x11_surface, Rectangle const& rectangle, cairo_content_t content, Color color, Window* window) :
  window_(window), rectangle_(rectangle)
{
  // Create an off-screen surface for tripple buffering.
  drawing_surface_ = cairo_surface_create_similar(x11_surface, content, rectangle.width(), rectangle.height());
  drawing_cr_ = cairo_create(drawing_surface_);

  if (content == CAIRO_CONTENT_COLOR)
  {
    // Set background color for off-screen surface.
    cairo_set_source_rgb(drawing_cr_, color.red(), color.green(), color.blue());

    // Fill surface with background color.
    cairo_paint(drawing_cr_);
  }
  else if (color.alpha() > 0)
  {
    // Set background color for off-screen surface.
    cairo_set_source_rgba(drawing_cr_, color.red(), color.green(), color.blue(), color.alpha());

    // Fill surface with background color.
    cairo_paint(drawing_cr_);
  }
}

Layer::~Layer()
{
  // Clean up.
  cairo_destroy(drawing_cr_);
  cairo_surface_destroy(drawing_surface_);
}

void Layer::draw(std::function<Rectangle(cairo_t*)> user_draw)
{
  cairo_save(drawing_cr_);
  // Apply layer offset, if any.
  cairo_translate(drawing_cr_, -offset_x(), -offset_y());

  // Call the drawing function of the user.
  Rectangle changed_area = user_draw(drawing_cr_);

  cairo_restore(drawing_cr_);

  window_->redraw(changed_area);
}

} // namespace cairowindow
