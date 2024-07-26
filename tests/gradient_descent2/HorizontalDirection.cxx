#include "sys.h"
#include "HorizontalDirection.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(HorizontalDirection hdirection)
{
  using enum HorizontalDirection;
  switch (hdirection)
  {
    AI_CASE_RETURN(left);
    AI_CASE_RETURN(undecided);
    AI_CASE_RETURN(right);
  }
  AI_NEVER_REACHED
}

std::string to_string(Region region)
{
  using enum Region;
  switch (region)
  {
    AI_CASE_RETURN(left);
    AI_CASE_RETURN(inbetween);
    AI_CASE_RETURN(right);
    AI_CASE_RETURN(unknown);
    AI_CASE_RETURN(invalid);
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
