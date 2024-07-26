#pragma once

#include <string>
#ifdef CWDEBUG
#include <ostream>
#endif

namespace gradient_descent {

// Is the next local extreme that we're looking for, a minimum or a maximum.
enum class ExtremeType
{
  minimum = -1,
  unknown = 0,
  maximum = 1
};

inline ExtremeType opposite(ExtremeType vdirection)
{
  return static_cast<ExtremeType>(-static_cast<int>(vdirection));
}

std::string to_string(ExtremeType vdirection);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, ExtremeType vdirection)
{
  return os << to_string(vdirection);
}
#endif

} // namespace gradient_descent
