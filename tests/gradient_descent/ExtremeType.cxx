#include "sys.h"
#include "ExtremeType.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(ExtremeType vdirection)
{
  switch (vdirection)
  {
    case ExtremeType::minimum: return "minimum";
    case ExtremeType::unknown: return "unknown";
    case ExtremeType::maximum: return "maximum";
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
