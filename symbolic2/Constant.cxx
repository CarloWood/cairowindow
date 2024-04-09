#include "sys.h"
#include "Constant.h"

namespace symbolic2 {

Expression const& Constant::s_cached_zero = get<Constant>(0);
Expression const& Constant::s_cached_one= get<Constant>(1);

bool Constant::equals(Expression const& other) const
{
  Constant const* other_constant = dynamic_cast<Constant const*>(&other);
  return other_constant && enumerator_ == other_constant->enumerator_ && denominator_ == other_constant->denominator_;
}

} // namespace symbolic2
