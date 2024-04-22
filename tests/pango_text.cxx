#include "cairowindow/FontPangoContext.h"
#include <cairo-svg.h>

#define WIDTH 1000
#define HEIGHT 500
#define TEXT_SIZE 13

int main(int argc, char* argv[])
{
  // Set the desired width and height of the SVG
  double width = WIDTH;
  double height = HEIGHT;

  char text[] = "The quick brown fox jumps over the lazy dog.";
  char output_filename[] = "image.png";

  //cr_surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, 1000, 300);
  cairo_surface_t* cr_surface = cairo_svg_surface_create("output.svg", width, height);

  cairo_t* cr = cairo_create(cr_surface);

  cairo_set_source_rgb(cr, 1.0, 1.0, 1.0);
  cairo_rectangle(cr, 0, 0, WIDTH, HEIGHT);
  cairo_fill(cr);

  cairo_set_source_rgb(cr, 0.0, 0.0, 0.0);

  cairowindow::FontPangoContext pango_context(cr);

  int ypos = 0;
  char font_desc[32];
  for (int i = 0; i < 30; i += 2)
  {
    sprintf(font_desc, "Ubuntu %u", i + 6);
    auto layout = pango_context.draw_text(font_desc, text);

    PangoRectangle pr;
    pango_layout_get_extents(layout, &pr, NULL);
    pango_extents_to_pixels(&pr, NULL);
    printf("x: %i, y: %i, width: %i, height: %i\n", pr.x, pr.y, pr.width, pr.height);

    ypos = ypos + pr.height + 2;
    if (ypos < HEIGHT)
      cairo_move_to(cr, 0, ypos);
    else
      break;
  }

  cairo_surface_write_to_png(cr_surface, output_filename);
  cairo_destroy(cr);
  cairo_surface_destroy(cr_surface);
}
