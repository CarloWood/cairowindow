#include "sys.h"
#include "Shape.h"
#ifdef CWDEBUG
#include "cairowindow/debug_channel.h"
#include "cairowindow/debugcairo.h"
#endif

namespace cairowindow::draw {

namespace {
int shape_index = 0;
} // namespace

ShapeEnum next_shape()
{
  return static_cast<ShapeEnum>(shape_index++ % number_of_shapes);
}

StrokeExtents Shape::do_draw(cairo_t* cr)
{
  DoutEntering(dc::cairowindow, "draw::Shape::do_draw(cr) [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
  Dout(dc::cairowindow, "geometry_ = " << geometry_ << "; style_ = " << style_);

  bool do_stroke = !style_.line_color().is_transparent();
  bool do_fill = !style_.fill_color().is_transparent();
  ASSERT(do_stroke || do_fill);

  // Draw the shape path.
  cairo_save(cr);
  if (style_.position() == at_corner)
  {
    // Center the shape around a corner of geometry [//].
    //
    //   offset_x
    //      |
    //      v
    // .----+----.
    // |    |    |
    // |    |    |
    // +----C----+ <-- offset_y
    // |    |////|  ^
    // |    |////|  |-- height
    // '----+----'  v             C: center of the shape
    //       <-->
    //       width
    cairo_translate(cr, geometry_.offset_x(), geometry_.offset_y());
    cairo_scale(cr, geometry_.width(), geometry_.height());
  }
  else if (style_.position() == at_center)
  {
    // Use geometry [//] as the border around the shape.
    //
    // offset_x
    // |
    // v
    // .----+----. <-- offset_y
    // |////|////|  ^
    // |////|////|  |
    // +----C----+  |-- height
    // |////|////|  |
    // |////|////|  |
    // '----+----'  v
    //  <------->
    //    width
    cairo_translate(cr, geometry_.offset_x() + 0.5 * geometry_.width(), geometry_.offset_y() + 0.5 * geometry_.height());
    cairo_scale(cr, 0.5 * geometry_.width(), 0.5 * geometry_.height());
  }
  else if (style_.position() == at_tip)
  {
    cairo_translate(cr, geometry_.offset_x() + geometry_.width(), geometry_.offset_y() + 0.5 * geometry_.height());
    cairo_rotate(cr, rotation_);
    cairo_translate(cr, -arrow_overshoot_, 0.0);
    cairo_scale(cr, geometry_.width(), geometry_.height());
  }
  switch (style_.shape())
  {
    case rectangle:
      //
      // -1,-1 ----> 1,-1
      //   ^          |
      //   |          |
      //   |          |
      //   |          v
      // -1,1 <----- 1,1
      cairo_move_to(cr, -1.0, -1.0);
      cairo_line_to(cr,  1.0, -1.0);
      cairo_line_to(cr,  1.0,  1.0);
      cairo_line_to(cr, -1.0,  1.0);
      cairo_close_path(cr);
      break;
    case ellipse:
      cairo_arc(cr, 0.0, 0.0, 1.0, 0.0, 2 * M_PI);
      break;
    case triangle_up:
      cairo_move_to(cr, -1.0,  1.0);
      cairo_line_to(cr,  0.0, -1.0);
      cairo_line_to(cr,  1.0,  1.0);
      cairo_close_path(cr);
      break;
    case triangle_down:
      cairo_move_to(cr, -1.0, -1.0);
      cairo_line_to(cr,  1.0, -1.0);
      cairo_line_to(cr,  0.0,  1.0);
      cairo_close_path(cr);
      break;
    case triangle_left:
      cairo_move_to(cr,  1.0, -1.0);
      cairo_line_to(cr,  1.0,  1.0);
      cairo_line_to(cr, -1.0,  0.0);
      cairo_close_path(cr);
      break;
    case triangle_right:
      cairo_move_to(cr, -1.0, -1.0);
      cairo_line_to(cr,  1.0,  0.0);
      cairo_line_to(cr, -1.0,  1.0);
      cairo_close_path(cr);
      break;
    case rhombus:
      cairo_move_to(cr,  0.0, -1.0);
      cairo_line_to(cr,  1.0,  0.0);
      cairo_line_to(cr,  0.0,  1.0);
      cairo_line_to(cr, -1.0,  0.0);
      cairo_close_path(cr);
      break;
    case diamond:
      cairo_move_to(cr, -0.66666, -1.0);
      cairo_line_to(cr,  0.66666, -1.0);
      cairo_line_to(cr,  1.0, -0.33333);
      cairo_line_to(cr,  0.0, 1.0);
      cairo_line_to(cr, -1.0, -0.33333);
      cairo_close_path(cr);
      break;
    case star:
    case cross:
      cairo_move_to(cr, -1.0, -1.0);
      cairo_line_to(cr,  1.0,  1.0);
      cairo_move_to(cr, -1.0,  1.0);
      cairo_line_to(cr,  1.0, -1.0);
      if (style_.shape() == cross)
        break;
      [[fallthrough]];
    case plus:
      cairo_move_to(cr,  0.0, -1.0);
      cairo_line_to(cr,  0.0,  1.0);
      cairo_move_to(cr, -1.0,  0.0);
      cairo_line_to(cr,  1.0,  0.0);
      break;
    case none_arrow_shape:
      // Why would we even get here?
      ASSERT(false);
      break;
    case open_arrow_shape:
      do_fill = false;
      [[fallthrough]];
    case filled_arrow_shape:
      // Draw arrows such that the tip is at 0,0.
      cairo_move_to(cr, -1.0, 0.5);
      cairo_line_to(cr, 0.0, 0.0);
      cairo_line_to(cr, -1.0, -0.5);
      if (do_fill)
        cairo_close_path(cr);
      break;
    case diamond_arrow_shape:
      cairo_move_to(cr, 0.0, 0.0);
      cairo_line_to(cr, -0.5, -0.5);
      cairo_line_to(cr, -1.0, 0.0);
      cairo_line_to(cr, -0.5, 0.5);
      cairo_close_path(cr);
      break;
    case circle_arrow_shape:
      cairo_arc(cr, -0.5, 0.0, 0.5, 0.0, 2 * M_PI);
      break;
  }
  cairo_restore(cr);

  double x1, y1;
  double x2, y2;

  if (do_fill)
    cairo_set_source_rgba(cr, style_.fill_color().red(), style_.fill_color().green(), style_.fill_color().blue(), style_.fill_color().alpha());
  if (do_stroke)
  {
    if (do_fill)
      cairo_fill_preserve(cr);
    cairo_set_source_rgba(cr, style_.line_color().red(), style_.line_color().green(), style_.line_color().blue(), style_.line_color().alpha());
    cairo_set_line_width(cr, style_.line_width());
    cairo_stroke_extents(cr, &x1, &y1, &x2, &y2);
    cairo_stroke(cr);
  }
  else
  {
    cairo_fill_extents(cr, &x1, &y1, &x2, &y2);
    cairo_fill(cr);
  }

  return {x1, y1, x2, y2};
}

} // namespace cairowindow::draw
