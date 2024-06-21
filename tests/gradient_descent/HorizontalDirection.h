#pragma once

#include <string>
#ifdef CWDEBUG
#include <ostream>
#endif

namespace gradient_descent {

enum class HorizontalDirection
{
  left = -1,
  undecided = 0,
  right = 1
};

inline HorizontalDirection opposite(HorizontalDirection hdirection)
{
  return static_cast<HorizontalDirection>(-static_cast<int>(hdirection));
}

std::string to_string(HorizontalDirection hdirection);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, HorizontalDirection hdirection)
{
  return os << to_string(hdirection);
}
#endif

} // namespace gradient_descent
