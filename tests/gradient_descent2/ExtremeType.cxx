#include "sys.h"
#include "ExtremeType.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(ExtremeType vdirection)
{
  using enum ExtremeType;
  switch (vdirection)
  {
    AI_CASE_RETURN(minimum);
    AI_CASE_RETURN(unknown);
    AI_CASE_RETURN(maximum);
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
