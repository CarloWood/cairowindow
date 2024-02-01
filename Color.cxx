#include "sys.h"
#include "Color.h"
#include "utils/ColorPool.h"

namespace cairowindow {

namespace {

utils::ColorPool<32> color_pool;

Color rgb(int r, int g, int b)
{
  return { r / 255.0, g / 255.0, b / 255.0 };
}

std::array<Color, 32> blackbg_color_table = {
  rgb(237, 21, 21),
  rgb(215, 255, 0),
  rgb(17, 209, 22),
  rgb(0, 135, 255),
  rgb(255, 95, 215),
  rgb(246, 116, 0),
  rgb(175, 95, 0),
  rgb(255, 255, 255),
  rgb(95, 255, 135),
  rgb(95, 95, 255),
  rgb(255, 0, 215),
  rgb(255, 95, 0),
  rgb(192, 57, 43),
  rgb(127, 140, 141),
  rgb(215, 0, 0),
  rgb(255, 215, 135),
  rgb(0, 215, 215),
  rgb(155, 89, 182),
  rgb(255, 135, 135),
  rgb(0, 135, 0),
  rgb(253, 188, 75),
  rgb(255, 215, 255),
  rgb(175, 0, 0),
  rgb(28, 220, 154),
  rgb(29, 153, 243),
  rgb(175, 135, 215),
  rgb(175, 135, 0),
  rgb(188, 188, 188),
  rgb(26, 188, 156),
  rgb(0, 95, 255),
  rgb(135, 255, 255),
  rgb(255, 215, 0)
};

std::array<Color, 32> whitebg_color_table = {
  rgb(97, 35, 2),
  rgb(212, 149, 30),
  rgb(0, 255, 0),
  rgb(72, 206, 255),
  rgb(0, 77, 215),
  rgb(255, 105, 180),
  rgb(248, 131, 121),
  rgb(16, 12, 8),
  rgb(34, 139, 34),
  rgb(0, 208, 201),
  rgb(0, 127, 255),
  rgb(255, 0, 255),
  rgb(231, 0, 0),
  rgb(211, 86, 0),
  rgb(59, 60, 54),
  rgb(30, 93, 31),
  rgb(54, 69, 79),
  rgb(85, 26, 139),
  rgb(204, 10, 87),
  rgb(226, 76, 0),
  rgb(237, 145, 33),
  rgb(173, 223, 173),
  rgb(0, 128, 128),
  rgb(0, 33, 71),
  rgb(169, 0, 163),
  rgb(185, 40, 28),
  rgb(128, 70, 27),
  rgb(141, 182, 0),
  rgb(1, 50, 32),
  rgb(102, 153, 204),
  rgb(178, 132, 190),
  rgb(128, 0, 0)
};

} // namespace

//static
Color Color::get_color(int color_index)
{
  return whitebg_color_table[color_index];
}

//static
Color Color::next_color()
{
  int color_index = color_pool.get_and_use_color();
  return whitebg_color_table[color_index];
}

} // namespace cairowindow
