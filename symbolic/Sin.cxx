#include "sys.h"
#include "Sin.h"
#include "Cos.h"
#include "Product.h"
#include <cmath>

namespace symbolic {

double Sin::evaluate() const
{
  return std::sin(arg_.evaluate());
}

Expression const& Sin::derivative(Symbol const& symbol) const
{
  return Product::multiply(arg_.derivative(symbol), Cos::realize(arg_));
}

} // namespace symbolic
