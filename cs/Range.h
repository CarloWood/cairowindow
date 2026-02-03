#pragma once

#include "CS.h"
#include "cairowindow/Range.h"
#ifdef CWDEBUG
#include "utils/to_string.h"
#endif

namespace cairowindow::cs {

template<CS cs>
class Range : public cairowindow::Range
{
 public:
  using cairowindow::Range::Range;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":";
    cairowindow::Range::print_on(os);
  }
#endif
};

} // namespace cairowindow::cs
