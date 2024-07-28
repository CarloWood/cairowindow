#include "sys.h"
#include "CubicToNextSampleType.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(CubicToNextSampleType state)
{
  using enum CubicToNextSampleType;
  switch (state)
  {
    AI_CASE_RETURN(unknown);
    AI_CASE_RETURN(flat);
    AI_CASE_RETURN(up);
    AI_CASE_RETURN(down);
    AI_CASE_RETURN(right_stop);
    AI_CASE_RETURN(left_stop);
    AI_CASE_RETURN(right_min);
    AI_CASE_RETURN(left_min);
    AI_CASE_RETURN(right_max);
    AI_CASE_RETURN(left_max);
    AI_CASE_RETURN(right_max_left_min);
    AI_CASE_RETURN(right_min_left_max);
    AI_CASE_RETURN(min);
    AI_CASE_RETURN(max);
    AI_CASE_RETURN(min_max);
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
