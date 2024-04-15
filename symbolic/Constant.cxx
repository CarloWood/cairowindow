#include "sys.h"
#include "Constant.h"

namespace symbolic {

bool Constant::equals(Expression const& other) const
{
  Constant const* other_constant = dynamic_cast<Constant const*>(&other);
  return other_constant && enumerator_ == other_constant->enumerator_ && denominator_ == other_constant->denominator_;
}

} // namespace symbolic
