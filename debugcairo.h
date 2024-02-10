#pragma once

#include <cairo/cairo.h>
#include "debug.h"

#ifdef CWDEBUG
NAMESPACE_DEBUG_CHANNELS_START
extern channel_ct cairo;
NAMESPACE_DEBUG_CHANNELS_END
#endif

namespace debugcairo {

std::ostream& operator<<(std::ostream& os, cairo_t const* cr);
std::ostream& operator<<(std::ostream& os, cairo_surface_t const* surface);
std::ostream& operator<<(std::ostream& os, cairo_text_extents_t const* text_extents);
std::ostream& operator<<(std::ostream& os, cairo_matrix_t const* matrix);

void debug_cairo_arc(cairo_t* cr, double xc, double yc, double radius, double angle1, double angle2);
void debug_cairo_clip(cairo_t* cr);
void debug_cairo_close_path(cairo_t* cr);
cairo_t* debug_cairo_create(cairo_surface_t* target COMMA_DEBUG_ONLY(std::string name));
void debug_cairo_destroy(cairo_t* cr);
void debug_cairo_fill(cairo_t* cr);
void debug_cairo_fill_extents(cairo_t* cr, double* x1, double* y1, double* x2, double* y2);
void debug_cairo_fill_preserve(cairo_t* cr);
void debug_cairo_get_matrix(cairo_t* cr, cairo_matrix_t* matrix);
cairo_surface_t* debug_cairo_get_target(cairo_t* cr);
void debug_cairo_line_to(cairo_t* cr, double x, double y);
void debug_cairo_matrix_transform_point(const cairo_matrix_t* matrix, double* x, double* y);
void debug_cairo_move_to(cairo_t* cr, double x, double y);
void debug_cairo_paint(cairo_t* cr);
cairo_pattern_t* debug_cairo_pattern_create_for_surface(cairo_surface_t* surface);
void debug_cairo_pattern_destroy(cairo_pattern_t* pattern);
void debug_cairo_pattern_set_extend(cairo_pattern_t* pattern, cairo_extend_t extend);
void debug_cairo_rectangle(cairo_t* cr, double x, double y, double width, double height);
void debug_cairo_restore(cairo_t* cr);
void debug_cairo_rotate(cairo_t* cr, double angle);
void debug_cairo_save(cairo_t* cr);
void debug_cairo_scale(cairo_t* cr, double sx, double sy);
void debug_cairo_select_font_face(cairo_t* cr, char const* family, cairo_font_slant_t slant, cairo_font_weight_t weight);
void debug_cairo_set_dash(cairo_t* cr, double const* dashes, int num_dashes, double offset);
void debug_cairo_set_font_size(cairo_t* cr, double size);
void debug_cairo_set_line_cap(cairo_t* cr, cairo_line_cap_t line_cap);
void debug_cairo_set_line_width(cairo_t* cr, double width);
void debug_cairo_set_operator(cairo_t* cr, cairo_operator_t op);
void debug_cairo_set_source(cairo_t* cr, cairo_pattern_t* source);
void debug_cairo_set_source_rgb(cairo_t* cr, double red, double green, double blue);
void debug_cairo_set_source_rgba(cairo_t* cr, double red, double green, double blue, double alpha);
void debug_cairo_set_source_surface(cairo_t* cr, cairo_surface_t* surface, double x, double y);
void debug_cairo_show_text(cairo_t* cr, char const* utf8);
void debug_cairo_stroke(cairo_t* cr);
void debug_cairo_stroke_extents(cairo_t* cr, double* x1, double* y1, double* x2, double* y2);
cairo_surface_t* debug_cairo_surface_create_similar(cairo_surface_t* other, cairo_content_t content, int width, int height
    COMMA_CWDEBUG_ONLY(std::string name));
void debug_cairo_surface_destroy(cairo_surface_t* surface);
void debug_cairo_text_extents(cairo_t* cr, char const* utf8, cairo_text_extents_t* extents);
void debug_cairo_translate(cairo_t* cr, double tx, double ty);
// Avoid using X11 macro "types", so we don't have to include the X11 headers.
cairo_surface_t* debug_cairo_xlib_surface_create(
    /*Display*/void* dpy,
    /*Drawable*/unsigned long d,
    /*Visual*/void* visual,
    int width, int height
    COMMA_CWDEBUG_ONLY(std::string name));

} // namespace debugcairo

#define cairo_arc(cr, xc, yc, radius, angle1, angle2) \
  debug_cairo_arc(cr, xc, yc, radius, angle1, angle2)

#define cairo_clip(cr) \
  debug_cairo_clip(cr)

#define cairo_close_path(cr) \
  debug_cairo_close_path(cr)

#define cairo_create(target) \
  debug_cairo_create(target)

#define cairo_destroy(cr) \
  debug_cairo_destroy(cr)

#define cairo_fill(cr) \
  debug_cairo_fill(cr)

#define cairo_fill_extents(cr, x1, y1, x2, y2) \
  debug_cairo_fill_extents(cr, x1, y1, x2, y2)

#define cairo_fill_preserve(cr) \
  debug_cairo_fill_preserve(cr)

#define cairo_get_matrix(cr, matrix) \
  debug_cairo_get_matrix(cr, matrix)

#define cairo_get_target(cr) \
  debug_cairo_get_target(cr)

#define cairo_line_to(cr, x, y) \
  debug_cairo_line_to(cr, x, y)

#define cairo_matrix_transform_point(matrix, x, y) \
  debug_cairo_matrix_transform_point(matrix, x, y)

#define cairo_move_to(cr, x, y) \
  debug_cairo_move_to(cr, x, y)

#define cairo_paint(cr) \
  debug_cairo_paint(cr)

#define cairo_pattern_create_for_surface(surface) \
  debug_cairo_pattern_create_for_surface(surface)

#define cairo_pattern_destroy(pattern) \
  debug_cairo_pattern_destroy(pattern)

#define cairo_pattern_set_extend(pattern, extend) \
  debug_cairo_pattern_set_extend(pattern, extend)

#define cairo_rectangle(cr, x, y, width, height) \
  debug_cairo_rectangle(cr, x, y, width, height)

#define cairo_restore(cr) \
  debug_cairo_restore(cr)

#define cairo_rotate(cr, angle) \
  debug_cairo_rotate(cr, angle)

#define cairo_save(cr) \
  debug_cairo_save(cr)

#define cairo_scale(cr, sx, sy) \
  debug_cairo_scale(cr, sx, sy)

#define cairo_select_font_face(cr, family, slant, weight) \
  debug_cairo_select_font_face(cr, family, slant, weight)

#define cairo_set_dash(cr, dashes, num_dashes, offset) \
  debug_cairo_set_dash(cr, dashes, num_dashes, offset)

#define cairo_set_font_size(cr, size) \
  debug_cairo_set_font_size(cr, size)

#define cairo_set_line_cap(cr, line_cap) \
  debug_cairo_set_line_cap(cr, line_cap)

#define cairo_set_line_width(cr, width) \
  debug_cairo_set_line_width(cr, width)

#define cairo_set_operator(cr, op) \
  debug_cairo_set_operator(cr, op)

#define cairo_set_source(cr, source) \
  debug_cairo_set_source(cr, source)

#define cairo_set_source_rgb(cr, red, green, blue) \
  debug_cairo_set_source_rgb(cr, red, green, blue)

#define cairo_set_source_rgba(cr, red, green, blue, alpha) \
  debug_cairo_set_source_rgba(cr, red, green, blue, alpha)

#define cairo_set_source_surface(cr, surface, x, y) \
  debug_cairo_set_source_surface(cr, surface, x, y)

#define cairo_show_text(cr, utf8) \
  debug_cairo_show_text(cr, utf8)

#define cairo_stroke(cr) \
  debug_cairo_stroke(cr)

#define cairo_stroke_extents(cr, x1, y1, x2, y2) \
  debug_cairo_stroke_extents(cr, x1, y1, x2, y2)

#define cairo_surface_create_similar(other, content, width, height) \
  debug_cairo_surface_create_similar(other, content, width, height)

#define cairo_surface_destroy(surface) \
  debug_cairo_surface_destroy(surface)

#define cairo_text_extents(cr, utf8, extents) \
  debug_cairo_text_extents(cr, utf8, extents)

#define cairo_translate(cr, tx, ty) \
  debug_cairo_translate(cr, tx, ty)

#define cairo_xlib_surface_create(dpy, d, visual, width, height) \
  debug_cairo_xlib_surface_create(dpy, d, visual, width, height)

namespace cairowindow {
using debugcairo::operator<<;
} // namespace cairowindow
