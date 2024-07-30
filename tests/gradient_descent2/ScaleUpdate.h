#pragma once

#include <string>
#ifdef CWDEBUG
#include <iostream>
#endif

namespace gradient_descent {

enum class ScaleUpdate
{
  first_sample,                 // When there is only one sample.
  initialized,                  // When there are two samples for the first time.
  towards_cp,                   // When we already had two samples (and therefore an "old" cubic approximation)
                                // and the new sample is at the critical point of the new cubic.
  away_from_cp,                 // When we already had two samples (and therefore an "old" cubic approximation)
                                // and the new sample was moved away from the critical point.
  disconnected                  // Same as away_from_vertex but the scale wasn't increased.
};

std::string to_string(ScaleUpdate scale_update);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, ScaleUpdate scale_update)
{
  return os << to_string(scale_update);
}
#endif

} // namespace gradient_descent
