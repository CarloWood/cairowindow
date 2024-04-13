#include "sys.h"
#include "Product.h"
#include "Power.h"
#include "debug.h"

namespace symbolic2 {

Power::Power(Expression const& arg1, Expression const& arg2) : BinaryOperator<Power>(arg1, arg2)
{
  // Can only raise expressions to constant powers.
  ASSERT(Constant::is_constant(arg2));
  // Don't use Power when arg2 is zero or one.
  ASSERT(!Constant::is_zero(arg2) && !Constant::is_one(arg2));
  // arg1 is not allowed to be a constant.
  ASSERT(!Constant::is_constant(arg1));
  // arg1 is not allowed to contain a constant factor.
  ASSERT(!Product::has_constant_factor(arg1));
}

//static
Expression const& Power::make_power(Expression const& base, Expression const& exponent)
{
  ASSERT(!Constant::is_constant(base));

  if (Constant::is_zero(exponent))
    return Constant::s_cached_one;
  else if (Constant::is_one(exponent))
    return base;

  return realize(base, exponent);
}

Expression const& Power::differentiate(Symbol const& symbol) const
{
  return Product::multiply(Product::multiply(get_exponent(), make_power(arg1_, get_exponent() + Constant::s_cached_minus_one)), arg1_.differentiate(symbol));
}

#ifdef SYMBOLIC2_PRINTING
void Power::print_on(std::ostream& os) const
{
  bool need_parens = needs_parens(arg1_.precedence(), Precedence::power, before);
  if (need_parens)
    os << '(';
  arg1_.print_on(os);
  if (need_parens)
    os << ')';
  os << '^';
  need_parens = needs_parens(arg2_.precedence(), Precedence::power, after);
  if (need_parens)
    os << '(';
  arg2_.print_on(os);
  if (need_parens)
    os << ')';
}
#endif

} // namespace symbolic2
