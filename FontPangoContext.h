#pragma once

#include <pango/pangocairo.h>
#include <string>

namespace cairowindow {

class FontPangoContext
{
 private:
  cairo_t* cr_;
  PangoContext* pc_;

  struct LayoutPtr
  {
    PangoLayout* ptr_;
    LayoutPtr(PangoContext* pc, std::string const& font_desc) : ptr_(pango_layout_new(pc))
    {
      PangoFontDescription* desc = pango_font_description_from_string(font_desc.c_str());
      pango_layout_set_font_description(ptr_, desc);
      pango_font_description_free(desc);
    }
    ~LayoutPtr() { g_object_unref(ptr_); }
    operator PangoLayout*() const { return ptr_; }
  };

 public:
  FontPangoContext(cairo_t* cr);

  LayoutPtr draw_text(std::string const& font_desc, std::string const& text);
};

} // namespace cairowindow
