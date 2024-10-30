#include "sys.h"
#include "Point.h"
#include "draw/Point.h"

namespace cairowindow {
namespace plot {

cairowindow::Rectangle const& Point::geometry() const
{
  return draw_object_->geometry();
}

} // namespace plot
} // namespace cairowindow
