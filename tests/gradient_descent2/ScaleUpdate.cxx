#include "sys.h"
#include "ScaleUpdate.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(ScaleUpdate scale_update)
{
  switch (scale_update)
  {
    AI_CASE_RETURN(ScaleUpdate::first_sample);
    AI_CASE_RETURN(ScaleUpdate::initialized);
    AI_CASE_RETURN(ScaleUpdate::towards_cp);
    AI_CASE_RETURN(ScaleUpdate::away_from_cp);
    AI_CASE_RETURN(ScaleUpdate::disconnected);
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
