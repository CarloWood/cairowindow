#include "sys.h"
#include "VerticalDirection.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(VerticalDirection vdirection)
{
  switch (vdirection)
  {
    case VerticalDirection::down: return "down";
    case VerticalDirection::unknown: return "unknown";
    case VerticalDirection::up: return "up";
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
