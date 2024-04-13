#include "sys.h"
#include "Product.h"
#include "Sum.h"
#include "Constant.h"
#include "Symbol.h"
#include "debug.h"

namespace symbolic2 {

Expression const& Product::differentiate(Symbol const& symbol) const
{
  return Sum::add(multiply(arg1_.differentiate(symbol), arg2_), multiply(arg1_, arg2_.differentiate(symbol)));
}

//static
Expression const& Product::make_product(Expression const& arg1, Expression const& arg2)
{
  // Only call make_product with a non-product arg1.
  ASSERT(!arg1.is_product());
  // Only call make_product with arg1 < arg2.
  ASSERT(is_less(arg1, arg2));

  if (Constant::is_zero(arg1))
    return arg1;
  else if (Constant::is_one(arg1))
    return arg2;

  return realize(arg1, arg2);
}

//static
Expression const& Product::combine(Expression const& arg1, Expression const& arg2)
{
  Expression const& base1 = arg1.get_base();
  Expression const& base2 = arg2.get_base();

  // Only call Product::combine for two powers of the same base.
  ASSERT(&base1 == &base2);

  return Power::make_power(base1, arg1.get_exponent() + arg2.get_exponent());
}

//static
Expression const& Product::multiply(Expression const& arg1, Expression const& arg2)
{
  ExpressionType type1 = arg1.type();
  ExpressionType type2 = arg2.type();

  if (type1 == sumT)
  {
    Sum const& sum1 = static_cast<Sum const&>(arg1);
    if (type2 == sumT)
    {
      Sum const& sum2 = static_cast<Sum const&>(arg2);
      return Sum::add(multiply(sum1.arg1(), sum2.arg1()), Sum::add(multiply(sum1.arg1(), sum2.arg2()),
            Sum::add(multiply(sum1.arg2(), sum2.arg1()), multiply(sum1.arg2(), sum2.arg2()))));
    }
    return Sum::add(multiply(sum1.arg1(), arg2), multiply(sum1.arg2(), arg2));
  }
  else if (type2 == sumT)
  {
    Sum const& sum2 = static_cast<Sum const&>(arg2);
    return Sum::add(multiply(arg1, sum2.arg1()), multiply(arg1, sum2.arg2()));
  }

  if (type1 == productT)
  {
    Product const* product1 = dynamic_cast<Product const*>(&arg1);
    if (type2 == productT)
    {
      Product const* product2 = dynamic_cast<Product const*>(&arg2);
      if (is_less(product1->arg1(), product2->arg1()))
        return make_product(product1->arg1(), multiply(product1->arg2(), arg2));
      else if (is_less(product2->arg1(), product1->arg1()))
        return make_product(product2->arg1(), multiply(arg1, product2->arg2()));
      else
        return multiply(multiply(combine(product1->arg1(), product2->arg1()), product1->arg2()), product2->arg2());
    }
    if (is_less(product1->arg1(), arg2))
      return make_product(product1->arg1(), multiply(product1->arg2(), arg2));
    else if (is_less(arg2, product1->arg1()))
      return make_product(arg2, arg1);
    else
      return multiply(combine(product1->arg1(), arg2), product1->arg2());
  }
  if (type2 == productT)
  {
    Product const* product2 = dynamic_cast<Product const*>(&arg2);
    if (is_less(arg1, product2->arg1()))
      return make_product(arg1, arg2);
    else if (is_less(product2->arg1(), arg1))
      return make_product(product2->arg1(), multiply(arg1, product2->arg2()));
    else
      return multiply(combine(arg1, product2->arg1()), product2->arg2());
  }

  // Neither is a Sum or Product.
  if (is_less(arg1, arg2))
    return make_product(arg1, arg2);
  else if (is_less(arg2, arg1))
    return make_product(arg2, arg1);
  else if (Constant const* constant1 = dynamic_cast<Constant const*>(&arg1))
  {
    ASSERT(Constant::is_constant(arg2));
    Constant const& constant2 = static_cast<Constant const&>(arg2);
    return Constant::realize(constant1->enumerator_ * constant2.enumerator_, constant1->denominator_ * constant2.denominator_);
  }
  else
    return combine(arg1, arg2);
}

//static
Expression const& Product::negate(Expression const& arg)
{
  Constant const& constant_factor = arg.get_constant_factor();
  if (Constant::is_one(constant_factor))
    return multiply(Constant::s_cached_minus_one, arg);

  // If arg has a constant factor that is not one, then it must be a product (get_constant_factor() called on a Constant return one).
  ASSERT(arg.is_product());
  Product const& product = static_cast<Product const&>(arg);
  return realize(Constant::s_cached_minus_one * static_cast<Constant const&>(product.arg1()), product.arg2());
}

#ifdef SYMBOLIC2_PRINTING
void Product::print_on(std::ostream& os) const
{
  Precedence my_precedence = precedence();
  if (my_precedence == Precedence::negation)
  {
    bool need_parens = needs_parens(arg2_.precedence(), Precedence::negation);
    os << "-";
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
  else
  {
    bool need_parens = needs_parens(arg1_.precedence(), my_precedence, before);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " * ";
    need_parens = needs_parens(arg2_.precedence(), my_precedence, after);
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
}
#endif

} // namespace symbolic2
