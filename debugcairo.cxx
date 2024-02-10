#include "sys.h"
#include "threadsafe/threadsafe.h"
#include "threadsafe/AIReadWriteSpinLock.h"
#include <cairo/cairo.h>
#include <cairo/cairo-xlib.h>
#include "debug.h"
#include <string>
#include <map>

#ifdef CWDEBUG
NAMESPACE_DEBUG_CHANNELS_START
channel_ct cairo("CAIRO");
NAMESPACE_DEBUG_CHANNELS_END
#endif

namespace debugcairo {

using cairo_t_pointer_map_container_t = std::map<cairo_t const*, std::string>;
using cairo_t_pointer_map_t = threadsafe::Unlocked<cairo_t_pointer_map_container_t, threadsafe::policy::ReadWrite<AIReadWriteSpinLock>>;
cairo_t_pointer_map_t cairo_t_pointer_map;

using cairo_surface_t_pointer_map_container_t = std::map<cairo_surface_t const*, std::string>;
using cairo_surface_t_pointer_map_t =
  threadsafe::Unlocked<cairo_surface_t_pointer_map_container_t, threadsafe::policy::ReadWrite<AIReadWriteSpinLock>>;
cairo_surface_t_pointer_map_t cairo_surface_t_pointer_map;

std::ostream& operator<<(std::ostream& os, cairo_t const* cr)
{
  cairo_t_pointer_map_t::rat cairo_t_pointer_map_r(cairo_t_pointer_map);
  auto it = cairo_t_pointer_map_r->find(cr);
  ASSERT(it != cairo_t_pointer_map_r->end());
  os << it->second;
  return os;
}

std::ostream& operator<<(std::ostream& os, cairo_surface_t const* surface)
{
  cairo_surface_t_pointer_map_t::rat cairo_surface_t_pointer_map_r(cairo_surface_t_pointer_map);
  auto it = cairo_surface_t_pointer_map_r->find(surface);
  ASSERT(it != cairo_surface_t_pointer_map_r->end());
  os << it->second;
  return os;
}

std::ostream& operator<<(std::ostream& os, cairo_text_extents_t const* text_extents)
{
  os << '{' <<
    "x_bearing:" << text_extents->x_bearing <<
    ", y_bearing:" << text_extents->y_bearing <<
    ", width:" << text_extents->width <<
    ", height:" << text_extents->height <<
    ", x_advance:" << text_extents->x_advance <<
    ", y_advance:" << text_extents->y_advance << '}';
  return os;
}

std::ostream& operator<<(std::ostream& os, cairo_matrix_t const* matrix)
{
  os <<
    "{xx:" << matrix->xx <<
    ", yx:" << matrix->yx <<
    ", xy:" << matrix->xy <<
    ", yy:" << matrix->yy <<
    ", x0:" << matrix->x0 <<
    ", y0:" << matrix->y0 << '}';
  return os;
}

void debug_cairo_arc(cairo_t* cr, double xc, double yc, double radius, double angle1, double angle2)
{
  Dout(dc::cairo, "cairo_arc(" << cr << ", " << xc << ", " << yc << ", " << radius << ", " << angle1 << ", " << angle2 << ")");
  cairo_arc(cr, xc, yc, radius, angle1, angle2);
}

void debug_cairo_clip(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_clip(" << cr << ")");
  cairo_clip(cr);
}

void debug_cairo_close_path(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_close_path(" << cr << ")");
  cairo_close_path(cr);
}

cairo_t* debug_cairo_create(cairo_surface_t* target COMMA_DEBUG_ONLY(std::string name))
{
  Dout(dc::cairo, "cairo_create(" << target << ")");
  cairo_t* cr = cairo_create(target);
  cairo_t_pointer_map_t::wat cairo_t_pointer_map_w(cairo_t_pointer_map);
  cairo_t_pointer_map_w->insert({cr, "cr:" + name});
  return cr;
}

void debug_cairo_destroy(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_destroy(" << cr << ")");
  cairo_destroy(cr);
}

void debug_cairo_fill(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_fill(" << cr << ")");
  cairo_fill(cr);
}

void debug_cairo_fill_extents(cairo_t* cr, double* x1, double* y1, double* x2, double* y2)
{
  Dout(dc::cairo|continued_cf, "cairo_fill_extents(" << cr << ", ");
  cairo_fill_extents(cr, x1, y1, x2, y2);
  Dout(dc::finish, "={" << *x1 << "}, ={" << *y1 << "}, ={" << *x2 << "}, ={" << *y2 << "})");
}

void debug_cairo_fill_preserve(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_fill_preserve(" << cr << ")");
  cairo_fill_preserve(cr);
}

void debug_cairo_get_matrix(cairo_t* cr, cairo_matrix_t* matrix)
{
  Dout(dc::cairo|continued_cf, "cairo_get_matrix(" << cr << ", ");
  cairo_get_matrix(cr, matrix);
  Dout(dc::finish, matrix << ")");
}

cairo_surface_t* debug_cairo_get_target(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_get_target(" << cr << ")");
  return cairo_get_target(cr);
}

void debug_cairo_line_to(cairo_t* cr, double x, double y)
{
  Dout(dc::cairo, "cairo_line_to(" << cr << ", " << x << ", " << y << ")");
  cairo_line_to(cr, x, y);
}

void debug_cairo_matrix_transform_point(cairo_matrix_t const* matrix, double* x, double* y)
{
  Dout(dc::cairo|continued_cf, "cairo_matrix_transform_point(" << matrix << ", ");
  cairo_matrix_transform_point(matrix, x, y);
  Dout(dc::finish, "={" << *x << "}, ={" << *y << "})");
}

void debug_cairo_move_to(cairo_t* cr, double x, double y)
{
  Dout(dc::cairo, "cairo_move_to(" << cr << ", " << x << ", " << y << ")");
  cairo_move_to(cr, x, y);
}

void debug_cairo_paint(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_paint(" << cr << ")");
  cairo_paint(cr);
}

cairo_pattern_t* debug_cairo_pattern_create_for_surface(cairo_surface_t* surface)
{
  Dout(dc::cairo, "cairo_pattern_create_for_surface(" << surface << ")");
  return cairo_pattern_create_for_surface(surface);
}

void debug_cairo_pattern_destroy(cairo_pattern_t* pattern)
{
  Dout(dc::cairo, "cairo_pattern_destroy(" << pattern << ")");
  cairo_pattern_destroy(pattern);
}

void debug_cairo_pattern_set_extend(cairo_pattern_t* pattern, cairo_extend_t extend)
{
  Dout(dc::cairo, "cairo_pattern_set_extend(" << pattern << ", " << extend << ")");
  cairo_pattern_set_extend(pattern, extend);
}

void debug_cairo_rectangle(cairo_t* cr, double x, double y, double width, double height)
{
  Dout(dc::cairo, "cairo_rectangle(" << cr << ", " << x << ", " << y << ", " << width << ", " << height << ")");
  ASSERT(width >= 0.0 && height >= 0.0);
  // If a negative value sneaks through, it is possible that it gets converted to
  // an unsigned short and ends up as 65536 - value.
  ASSERT(x >= 0.0 && y >= 0.0 && x < 32768.0 && y < 32768.0);
  cairo_rectangle(cr, x, y, width, height);
}

void debug_cairo_restore(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_restore(" << cr << ")");
  cairo_restore(cr);
}

void debug_cairo_rotate(cairo_t* cr, double angle)
{
  Dout(dc::cairo, "cairo_rotate(" << cr << ", " << angle << ")");
  cairo_rotate(cr, angle);
}

void debug_cairo_save(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_save(" << cr << ")");
  cairo_save(cr);
}

void debug_cairo_scale(cairo_t* cr, double sx, double sy)
{
  Dout(dc::cairo, "cairo_scale(" << cr << ", " << sx << ", " << sy << ")");
  cairo_scale(cr, sx, sy);
}

void debug_cairo_select_font_face(cairo_t* cr, char const* family, cairo_font_slant_t slant, cairo_font_weight_t weight)
{
  Dout(dc::cairo, "cairo_select_font_face(" << cr << ", \"" << family << "\", " << slant << ", " << weight << ")");
  cairo_select_font_face(cr, family, slant, weight);
}

void debug_cairo_set_dash(cairo_t* cr, double const* dashes, int num_dashes, double offset)
{
  Dout(dc::cairo, "cairo_set_dash(" << cr << ", " << dashes << ", " << num_dashes << ", " << offset << ")");
  cairo_set_dash(cr, dashes, num_dashes, offset);
}

void debug_cairo_set_font_size(cairo_t* cr, double size)
{
  Dout(dc::cairo, "cairo_set_font_size(" << cr << ", " << size << ")");
  cairo_set_font_size(cr, size);
}

void debug_cairo_set_line_cap(cairo_t* cr, cairo_line_cap_t line_cap)
{
  Dout(dc::cairo, "cairo_set_line_cap(" << cr << ", " << line_cap << ")");
  cairo_set_line_cap(cr, line_cap);
}

void debug_cairo_set_line_width(cairo_t* cr, double width)
{
  Dout(dc::cairo, "cairo_set_line_width(" << cr << ", " << width << ")");
  cairo_set_line_width(cr, width);
}

void debug_cairo_set_operator(cairo_t* cr, cairo_operator_t op)
{
  Dout(dc::cairo, "cairo_set_operator(" << cr << ", " << op << ")");
  cairo_set_operator(cr, op);
}

void debug_cairo_set_source(cairo_t* cr, cairo_pattern_t* source)
{
  Dout(dc::cairo, "cairo_set_source(" << cr << ", " << source << ")");
  cairo_set_source(cr, source);
}

void debug_cairo_set_source_rgb(cairo_t* cr, double red, double green, double blue)
{
  Dout(dc::cairo, "cairo_set_source_rgb(" << cr << ", " << red << ", " << green << ", " << blue << ")");
  cairo_set_source_rgb(cr, red, green, blue);
}

void debug_cairo_set_source_rgba(cairo_t* cr, double red, double green, double blue, double alpha)
{
  Dout(dc::cairo, "cairo_set_source_rgba(" << cr << ", " << red << ", " << green << ", " << blue << ", " << alpha << ")");
  cairo_set_source_rgba(cr, red, green, blue, alpha);
}

void debug_cairo_set_source_surface(cairo_t* cr, cairo_surface_t* surface, double x, double y)
{
  Dout(dc::cairo, "cairo_set_source_surface(" << cr << ", " << surface << ", " << x << ", " << y << ")");
  cairo_set_source_surface(cr, surface, x, y);
}

void debug_cairo_show_text(cairo_t* cr, char const* utf8)
{
  Dout(dc::cairo, "cairo_show_text(" << cr << ", \"" << utf8 << "\")");
  cairo_show_text(cr, utf8);
}

void debug_cairo_stroke(cairo_t* cr)
{
  Dout(dc::cairo, "cairo_stroke(" << cr << ")");
  cairo_stroke(cr);
}

void debug_cairo_stroke_extents(cairo_t* cr, double* x1, double* y1, double* x2, double* y2)
{
  Dout(dc::cairo|continued_cf, "cairo_stroke_extents(" << cr << ", ");
  cairo_stroke_extents(cr, x1, y1, x2, y2);
  Dout(dc::finish, "={" << *x1 << "}, ={" << *y1 << "}, ={" << *x2 << "}, ={" << *y2 << "})");
}

cairo_surface_t* debug_cairo_surface_create_similar(cairo_surface_t* other, cairo_content_t content, int width, int height
    COMMA_CWDEBUG_ONLY(std::string debug_name))
{
  Dout(dc::cairo, "cairo_surface_create_similar(" << other << ", " << content << ", " << width << ", " << height << ")");
  cairo_surface_t* surface = cairo_surface_create_similar(other, content, width, height);
  cairo_surface_t_pointer_map_t::wat cairo_surface_t_pointer_map_w(cairo_surface_t_pointer_map);
  cairo_surface_t_pointer_map_w->insert({surface, "surface:\"" + debug_name + '"'});
  return surface;
}

void debug_cairo_surface_destroy(cairo_surface_t* surface)
{
  Dout(dc::cairo, "cairo_surface_destroy(" << surface << ")");
  cairo_surface_destroy(surface);
}

void debug_cairo_text_extents(cairo_t* cr, char const* utf8, cairo_text_extents_t* extents)
{
  Dout(dc::cairo|continued_cf, "cairo_text_extents(" << cr << ", \"" << utf8 << "\", ");
  cairo_text_extents(cr, utf8, extents);
  Dout(dc::finish, "=" << extents << ")");
}

void debug_cairo_translate(cairo_t* cr, double tx, double ty)
{
  Dout(dc::cairo, "cairo_translate(" << cr << ", " << tx << ", " << ty << ")");
  cairo_translate(cr, tx, ty);
}

cairo_surface_t* debug_cairo_xlib_surface_create(
    /*Display*/void* dpy,
    /*Drawable*/unsigned long d,
    /*Visual*/void* visual,
    int width, int height
    COMMA_CWDEBUG_ONLY(std::string name))
{
  Display* x11_dpy = (Display*)dpy;
  Drawable x11_d = d;
  Visual* x11_visual = (Visual*)visual;
  Dout(dc::cairo, "cairo_xlib_surface_create(dpy, visual, " << width << ", " << height << ")");
  cairo_surface_t* surface = cairo_xlib_surface_create(x11_dpy, x11_d, x11_visual, width, height);
  cairo_surface_t_pointer_map_t::wat cairo_surface_t_pointer_map_w(cairo_surface_t_pointer_map);
  cairo_surface_t_pointer_map_w->insert({surface, "surface:" + name});
  return surface;
}

} // namespace debugcairo
