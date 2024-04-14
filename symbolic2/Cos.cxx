#include "sys.h"
#include "Cos.h"
#include "Sin.h"
#include "Product.h"
#include <cmath>

namespace symbolic2 {

double Cos::evaluate() const
{
  return std::cos(arg_.evaluate());
}

Expression const& Cos::differentiate(Symbol const& symbol) const
{
  return Product::multiply(Product::negate(arg_.differentiate(symbol)), Sin::realize(arg_));
}

} // namespace symbolic2
