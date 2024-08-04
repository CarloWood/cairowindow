#pragma once

#include <string>
#ifdef CWDEBUG
#include <iostream>
#endif

namespace gradient_descent {

enum class IterationState
{
  initialization,
  first_cubic,
  find_extreme,
  next_sample,
  finish,
  success
};

std::string to_string(IterationState state);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, IterationState state)
{
  return os << to_string(state);
}
#endif

} // namespace gradient_descent
