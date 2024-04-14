#include "sys.h"
#include "Sin.h"
#include "Cos.h"
#include "Product.h"
#include <cmath>

namespace symbolic2 {

double Sin::evaluate() const
{
  return std::sin(arg_.evaluate());
}

Expression const& Sin::differentiate(Symbol const& symbol) const
{
  return Product::multiply(arg_.differentiate(symbol), Cos::realize(arg_));
}

} // namespace symbolic2
