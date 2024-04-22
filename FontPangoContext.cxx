#include "sys.h"
#include "FontPangoContext.h"

namespace cairowindow {

FontPangoContext::FontPangoContext(cairo_t* cr) : cr_(cr)
{
  pc_ = pango_cairo_create_context(cr);
  cairo_font_options_t* font_options = cairo_font_options_create();

  cairo_font_options_set_antialias(font_options, CAIRO_ANTIALIAS_SUBPIXEL);
  cairo_font_options_set_subpixel_order(font_options, CAIRO_SUBPIXEL_ORDER_RGB);
  cairo_font_options_set_hint_style(font_options, CAIRO_HINT_STYLE_SLIGHT);

  pango_cairo_context_set_font_options(pc_, font_options);
  cairo_font_options_destroy(font_options);    // Free font description after use.
}

FontPangoContext::LayoutPtr FontPangoContext::draw_text(std::string const& font_desc, std::string const& text)
{
  LayoutPtr layout{pc_, font_desc};
  pango_layout_set_text(layout, text.c_str(), -1);
  pango_cairo_update_layout(cr_, layout);
  pango_cairo_show_layout(cr_, layout);
  return layout;
}

} // namespace cairowindow
