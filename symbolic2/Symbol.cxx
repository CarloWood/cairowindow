#include "sys.h"
#include "Constant.h"
#include "Symbol.h"

namespace symbolic2 {

bool Symbol::equals(Expression const& other) const
{
  Symbol const* other_symbol = dynamic_cast<Symbol const*>(&other);
  return other_symbol && name_ == other_symbol->name_;
}

} // namespace symbolic2
