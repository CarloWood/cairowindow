#include "sys.h"
#include "HorizontalDirection.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(HorizontalDirection2 hdirection)
{
  switch (hdirection)
  {
    case HorizontalDirection2::left: return "left";
    case HorizontalDirection2::undecided: return "undecided";
    case HorizontalDirection2::right: return "right";
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
  }
  AI_NEVER_REACHED
}

std::string to_string(Restriction restriction)
{
  switch (restriction)
  {
    case Restriction::left: return "left";
    case Restriction::none: return "none";
    case Restriction::right: return "right";
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
