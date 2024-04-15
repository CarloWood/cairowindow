#include "sys.h"
#include "Log.h"
#include "Product.h"
#include "Power.h"
#include <cmath>

namespace symbolic2 {

double Log::evaluate() const
{
  return std::log(arg_.evaluate());
}

Expression const& Log::derivative(Symbol const& symbol) const
{
  return Product::multiply(arg_.derivative(symbol), Power::make_power(arg_, Constant::s_cached_minus_one));
}

} // namespace symbolic2
