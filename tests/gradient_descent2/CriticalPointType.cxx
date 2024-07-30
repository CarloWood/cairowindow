#include "sys.h"
#include "CriticalPointType.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(CriticalPointType critical_point_type)
{
  switch (critical_point_type)
  {
    AI_CASE_RETURN(CriticalPointType::none);
    AI_CASE_RETURN(CriticalPointType::minimum);
    AI_CASE_RETURN(CriticalPointType::maximum);
    AI_CASE_RETURN(CriticalPointType::inflection_point);
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
