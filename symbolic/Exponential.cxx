#include "sys.h"
#include "Exponential.h"
#include "Product.h"
#include <cmath>

namespace symbolic {

double Exponential::evaluate() const
{
  return std::exp(arg_.evaluate());
}

Expression const& Exponential::derivative(Symbol const& symbol) const
{
  return Product::multiply(arg_.derivative(symbol), *this);
}

} // namespace symbolic
