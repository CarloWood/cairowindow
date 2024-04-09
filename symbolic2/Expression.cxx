#include "sys.h"
#include "Expression.h"
#include "debug.h"
#ifdef CWDEBUG
#include <iomanip>
#endif

namespace symbolic2 {

//static
std::unordered_set<std::unique_ptr<Expression>, std::hash<std::unique_ptr<Expression>>, KeyEqual> Expression::s_database_;

#ifdef CWDEBUG
//static
void Expression::dump_database()
{
  for (auto&& ptr : s_database_)
    Dout(dc::notice, std::hex << std::hash<std::unique_ptr<Expression>>{}(ptr) << " @ " << (void*)ptr.get() << " : " << *ptr);
}
#endif

} // namespace symbolic2
