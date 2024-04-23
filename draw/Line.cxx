#include "sys.h"
#include "Line.h"
#include "utils/macros.h"
#include <iostream>
#include "debug.h"

namespace cairowindow::draw {

#ifdef CWDEBUG
char const* to_string(LineCap line_cap)
{
  switch (line_cap)
  {
    AI_CASE_RETURN(LineCap::undefined);
    AI_CASE_RETURN(LineCap::butt);
    AI_CASE_RETURN(LineCap::round);
    AI_CASE_RETURN(LineCap::square);
  }
  AI_NEVER_REACHED
}

std::ostream& operator<<(std::ostream& os, LineCap line_cap)
{
  return os << to_string(line_cap);
}

#endif

} // namespace cairowindow::draw
