#include "sys.h"
#include "Product.h"
#include "Power.h"
#include "debug.h"

namespace symbolic {

Power::Power(Expression const& arg1, Expression const& arg2) : BinaryOperator<Power>(arg1, arg2)
{
  // Can only raise expressions to constant powers.
  ASSERT(arg2.is_constant());
  // Don't use Power when arg2 is zero or one.
  ASSERT(!Constant::is_zero(arg2) && !Constant::is_one(arg2));
  // arg1 is not allowed to be a constant.
  ASSERT(!arg1.is_constant());
  // arg1 is not allowed to be a product.
  ASSERT(!arg1.is_product());
  // arg1 is not allowed to be a power.
  ASSERT(!arg1.is_power());
}

//static
Expression const& Power::make_power2(Expression const& base, Constant const& exponent)
{
  DoutEntering(dc::symbolic|continued_cf, "Power::make_power2(" << base << ", " << exponent << ") --> ");

  if (base.is_product())
  {
    Product const& product = static_cast<Product const&>(base);
    if (product.arg1().is_power())
    {
      Power const& power = static_cast<Power const&>(product.arg1());
      Expression const& result = Product::multiply(make_power(power.get_base(), power.get_exponent() * exponent), make_power2(product.arg2(), exponent));
      Dout(dc::finish, result);
      return result;
    }
    if (product.arg1().is_constant())
    {
      Expression const& result = Product::multiply(product.arg1() ^ exponent, make_power2(product.arg2(), exponent));
      Dout(dc::finish, result);
      return result;
    }
    Expression const& result = Product::multiply(realize(product.arg1(), exponent), make_power2(product.arg2(), exponent));
    Dout(dc::finish, result);
    return result;
  }

  if (base.is_power())
  {
    Power const& power = static_cast<Power const&>(base);
    Expression const& result = make_power(power.get_base(), power.get_exponent() * exponent);
    Dout(dc::finish, result);
    return result;
  }

  Expression const& result = realize(base, exponent);
  Dout(dc::finish, result);
  return result;
}

//static
Expression const& Power::make_power(Expression const& base, Constant const& exponent)
{
  DoutEntering(dc::symbolic|continued_cf, "Power::make_power(" << base << ", " << exponent << ") --> ");

  ASSERT(!base.is_constant());

  if (Constant::is_zero(exponent))
  {
    Dout(dc::finish, "-1");
    return Constant::s_cached_one;
  }
  else if (Constant::is_one(exponent))
  {
    Dout(dc::finish, base);
    return base;
  }

  Expression const& result = make_power2(base, exponent);
  Dout(dc::finish, result);
  return result;
}

Expression const& Power::derivative(Symbol const& symbol) const
{
  return Product::multiply(Product::multiply(get_exponent(), make_power(arg1_, get_exponent() + Constant::s_cached_minus_one)), arg1_.derivative(symbol));
}

#ifdef SYMBOLIC_PRINTING
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

} // namespace symbolic
