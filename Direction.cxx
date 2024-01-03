#include "sys.h"
#include "Direction.h"
#include "Line.h"

namespace cairowindow {

Direction::Direction(Line const& line) : Direction(line.direction())
{
}

//static
Direction const Direction::up{0.0, 1.0};
//static
Direction const Direction::down{0.0, -1.0};
//static
Direction const Direction::left{-1.0, 0.0};
//static
Direction const Direction::right{1.0, 0.0};

} // namespace cairowindow
