#include "sys.h"
#include "Shape.h"

namespace cairowindow::draw {

namespace {
int shape_index = 0;
} // namespace

ShapeEnum next_shape()
{
  return static_cast<ShapeEnum>(shape_index++ % number_of_shapes);
}

#ifdef CWDEBUG
void ShapeStyle::print_on(std::ostream& os) const
{
  os << "{line_color:" << line_color <<
       ", fill_color:" << fill_color <<
       ", line_width:" << line_width <<
       ", at_corner:" << std::boolalpha << at_corner <<
       ", shape:" << shape << "}";
}
#endif

} // namespace cairowindow::draw
