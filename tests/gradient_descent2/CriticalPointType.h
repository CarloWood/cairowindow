#pragma once

#include <string>
#ifdef CWDEBUG
#include <iostream>
#endif

namespace gradient_descent {

enum class CriticalPointType
{
  none,
  minimum,
  maximum,
  inflection_point
};

std::string to_string(CriticalPointType critical_point_type);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, CriticalPointType critical_point_type)
{
  return os << to_string(critical_point_type);
}
#endif

} // namespace gradient_descent
