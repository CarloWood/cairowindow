#include "sys.h"
#include "Direction.h"
#include "Line.h"

namespace cairowindow {

Direction::Direction(Line const& line) : Direction(line.normal().y(), -line.normal().x())
{
}

} // namespace cairowindow
