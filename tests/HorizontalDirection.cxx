#include "sys.h"
#include "HorizontalDirection.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(HorizontalDirection hdirection)
{
  switch (hdirection)
  {
    case HorizontalDirection::left: return "left";
    case HorizontalDirection::undecided: return "undecided";
    case HorizontalDirection::right: return "right";
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
