#include "sys.h"
#include "ArrowHead.h"

namespace cairowindow::draw {

//static
std::array<ArrowHead::Size, number_of_arrow_shapes> ArrowHead::s_arrow_head_size = {{
  {.width = 0, .height = 0},    // none_arrow
  {.width = 10, .height = 5},   // open_arrow
  {.width = 10, .height = 5},   // filled_arrow
  {.width = 20, .height = 5},   // diamond_arrow
  {.width = 10, .height =10}    // circle_arrow
}};

} // namespace cairowindow::draw
