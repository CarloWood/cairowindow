#include "sys.h"
#include "IterationState.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(IterationState state)
{
  using enum IterationState;
  switch (state)
  {
    AI_CASE_RETURN(initialization);
    AI_CASE_RETURN(first_cubic);
    AI_CASE_RETURN(find_extreme);
    AI_CASE_RETURN(next_sample);
    AI_CASE_RETURN(success);
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
