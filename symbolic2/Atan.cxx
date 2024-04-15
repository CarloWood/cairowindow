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

Expression const& Atan::derivative(Symbol const& symbol) const
{
  return Product::multiply(arg_.derivative(symbol), Constant::s_cached_one / Sum::add(Constant::s_cached_one, arg_ ^ 2));
}

} // namespace symbolic2
