#include "sys.h"
#include "Constant.h"
#include "Expression.h"
#include "debug.h"
#ifdef CWDEBUG
#include <iomanip>
#endif

namespace symbolic2 {

//static
std::unordered_set<std::unique_ptr<Expression>, std::hash<std::unique_ptr<Expression>>, KeyEqual> Expression::s_database_;

//static - these must be initialized AFTER Expression::s_database_.
Constant const& Constant::s_cached_zero = Constant::realize(0);
Constant const& Constant::s_cached_one = Constant::realize(1);
Constant const& Constant::s_cached_minus_one = Constant::realize(-1);

Constant const& Expression::get_constant_factor() const
{
  return Constant::s_cached_one;
}

Constant const& Expression::get_exponent() const
{
  return Constant::s_cached_one;
}

#ifdef CWDEBUG
//static
void Expression::dump_database()
{
  for (auto&& ptr : s_database_)
    Dout(dc::notice, std::hex << std::hash<std::unique_ptr<Expression>>{}(ptr) << " @ " << (void*)ptr.get() << " : " << std::dec << *ptr);
}
#endif

} // namespace symbolic2
