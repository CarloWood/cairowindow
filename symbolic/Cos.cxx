#include "sys.h"
#include "Cos.h"
#include "Sin.h"
#include "Product.h"
#include <cmath>

namespace symbolic {

double Cos::evaluate() const
{
  return std::cos(arg_.evaluate());
}

Expression const& Cos::derivative(Symbol const& symbol) const
{
  return Product::multiply(Product::negate(arg_.derivative(symbol)), Sin::realize(arg_));
}

} // namespace symbolic
