#include "sys.h"
#include "Atan.h"
#include "Sum.h"
#include "Product.h"
#include <cmath>

namespace symbolic2 {

double Atan::evaluate() const
{
  return std::atan(arg_.evaluate());
}

Expression const& Atan::differentiate(Symbol const& symbol) const
{
  return Product::multiply(arg_.differentiate(symbol), Constant::s_cached_one / Sum::add(Constant::s_cached_one, arg_ ^ 2));
}

} // namespace symbolic2
