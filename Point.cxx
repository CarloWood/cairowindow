#include "sys.h"
#include "Point.h"
#include "Direction.h"

namespace cairowindow {

Point Point::operator+(Direction const& direction)
{
  return {x_ + direction.x(), y_ + direction.y()};
}

} // namespace cairowindow
