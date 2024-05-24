#pragma once

#include <string>
#ifdef CWDEBUG
#include <ostream>
#endif

namespace gradient_descent {

enum class VerticalDirection
{
  down = -1,
  up = 1
};

inline VerticalDirection opposite(VerticalDirection vdirection)
{
  return static_cast<VerticalDirection>(-static_cast<int>(vdirection));
}

std::string to_string(VerticalDirection vdirection);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, VerticalDirection vdirection)
{
  return os << to_string(vdirection);
}
#endif

} // namespace gradient_descent
