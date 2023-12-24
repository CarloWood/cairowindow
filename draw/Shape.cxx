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

} // namespace cairowindow::draw
