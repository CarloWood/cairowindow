#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/Color.h"

namespace cairowindow::draw {

enum ShapeEnum
{
  rectangle,
  ellipse,
  triangle,
  triangle_up = triangle,
  triangle_down,
  triangle_left,
  triangle_right,
  rhombus,
  diamond,
  cross,
  plus,
  star
};

static constexpr int number_of_shapes = star + 1;

ShapeEnum next_shape();

struct ShapeStyle
{
  // Must start with the same style variables as CircleStyle, because we do a reinterpret_cast between the two!
  Color line_color = color::transparent;
  Color fill_color = color::transparent;
  double line_width = 2.0;
  bool at_corner = false;
  // This only exists in ShapeStyle.
  ShapeEnum shape = rectangle;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

class Shape : public LayerRegion
{
 protected:
  Rectangle geometry_;
  ShapeStyle style_;

 public:
  Shape(Rectangle const& geometry, ShapeStyle style) : geometry_(geometry), style_(style) { }

  ShapeStyle& style() { return style_; }
  ShapeStyle const& style() const { return style_; }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::notice, "draw::Shape::do_draw(cr) [" << this << "]");
    Dout(dc::notice, "geometry_ = " << geometry_ << "; style_ = " << style_);

    bool do_stroke = !style_.line_color.is_transparent();
    bool do_fill = !style_.fill_color.is_transparent();
    ASSERT(do_stroke || do_fill);

    // Draw the shape path.
    cairo_save(cr);
    if (style_.at_corner)
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
    else
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
    switch (style_.shape)
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
        if (style_.shape == cross)
          break;
        [[fallthrough]];
      case plus:
        cairo_move_to(cr,  0.0, -1.0);
        cairo_line_to(cr,  0.0,  1.0);
        cairo_move_to(cr, -1.0,  0.0);
        cairo_line_to(cr,  1.0,  0.0);
        break;
    }
    cairo_restore(cr);

    double x1, y1;
    double x2, y2;

    if (do_fill)
      cairo_set_source_rgba(cr, style_.fill_color.red(), style_.fill_color.green(), style_.fill_color.blue(), style_.fill_color.alpha());
    if (do_stroke)
    {
      if (do_fill)
        cairo_fill_preserve(cr);
      cairo_set_source_rgba(cr, style_.line_color.red(), style_.line_color.green(), style_.line_color.blue(), style_.line_color.alpha());
      cairo_set_line_width(cr, style_.line_width);
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
};

} // namespace cairowindow::draw
