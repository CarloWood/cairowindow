#include "sys.h"
#include "ChessPiece.h"

namespace cairowindow::draw {

StrokeExtents ChessPiece::do_draw(cairo_t* cr)
{
  DoutEntering(dc::cairowindow, "draw::ChessPiece::do_draw(cr) [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  //cairo_set_source_rgba(cr, style_.line_color().red(), style_.line_color().green(), style_.line_color().blue(), style_.line_color().alpha());
  cairo_set_source_rgba(cr, 1, 1, 1, 1);
//    cairo_set_line_width(cr, 2.0);
  cairo_rectangle(cr, x1_ + 10, y1_ + 10, square_size_ - 20, square_size_ - 20);
  double x1, y1, x2, y2;
  cairo_fill_extents(cr, &x1, &y1, &x2, &y2);
  cairo_fill(cr);

  return {x1, y1, x2, y2};
}

} // namespace cairowindow::draw
