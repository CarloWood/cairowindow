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

std::string to_string(Region region)
{
  switch (region)
  {
    case Region::left: return "left";
    case Region::inbetween: return "inbetween";
    case Region::right: return "right";
    case Region::unknown: return "unknown";
    case Region::invalid: return "invalid";
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
