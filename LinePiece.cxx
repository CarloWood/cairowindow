#include "sys.h"
#include "LinePiece.h"
#include "Direction.h"

namespace cairowindow {

Direction LinePiece::direction() const
{
  return Direction{from_, to_};
}

} // namespace cairowindow
