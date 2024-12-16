#include "sys.h"
#include "Product.h"
#include "Sum.h"
#include "Constant.h"
#include "Symbol.h"
#include "debug.h"

namespace symbolic {

Expression const& Product::derivative(Symbol const& symbol) const
{
  return Sum::add(multiply(arg1_.derivative(symbol), arg2_), multiply(arg1_, arg2_.derivative(symbol)));
}

//static
Expression const& Product::make_product(Expression const& arg1, Expression const& arg2)
{
  DoutEntering(dc::symbolic|continued_cf, "Product::make_product(" << arg1 << ", " << arg2 << ") --> ");

  // Only call make_product with a non-product arg1.
  ASSERT(!arg1.is_product());
  // Only call make_product with arg1 < arg2.
  ASSERT(is_less(arg1, arg2));

  if (Constant::is_zero(arg1))
  {
    Dout(dc::finish, arg1);
    return Constant::s_cached_zero;
  }
  else if (Constant::is_one(arg1))
  {
    Dout(dc::finish, arg2);
    return arg2;
  }

  Expression const& result = realize(arg1, arg2);
  Dout(dc::finish, result);
  return result;
}

//static
Expression const& Product::combine(Expression const& arg1, Expression const& arg2)
{
  DoutEntering(dc::symbolic|continued_cf, "Product::combine(" << arg1 << ", " << arg2 << ") --> ");

  Constant const* constant1 = dynamic_cast<Constant const*>(&arg1);
  if (constant1)
  {
    if (arg2.is_zero_function())
      return Constant::s_cached_zero;
    ASSERT(arg2.type() == constantT);
    Expression const& result = *constant1 * static_cast<Constant const&>(arg2);
    Dout(dc::finish, result);
    return result;
  }
  else if (arg1.is_zero_function())
    return Constant::s_cached_zero;

  Expression const& base1 = arg1.get_base();
  Expression const& base2 = arg2.get_base();

  // Only call Product::combine for two powers of the same base.
  ASSERT(&base1 == &base2);

  Expression const& result = Power::make_power(base1, arg1.get_exponent() + arg2.get_exponent());
  Dout(dc::finish, result);
  return result;
}

//static
Expression const& Product::multiply(Expression const& arg1, Expression const& arg2)
{
  DoutEntering(dc::symbolic|continued_cf, "Product::multiply(" << arg1 << ", " << arg2 << ") --> ");

  ExpressionType type1 = arg1.type();
  ExpressionType type2 = arg2.type();

  if (type1 == sumT)
  {
    Sum const& sum1 = static_cast<Sum const&>(arg1);
    if (type2 == sumT)
    {
      Sum const& sum2 = static_cast<Sum const&>(arg2);
      Expression const& result = Sum::add(multiply(sum1.arg1(), sum2.arg1()), Sum::add(multiply(sum1.arg1(), sum2.arg2()),
            Sum::add(multiply(sum1.arg2(), sum2.arg1()), multiply(sum1.arg2(), sum2.arg2()))));
      Dout(dc::finish, result);
      return result;
    }
    Expression const& result = Sum::add(multiply(sum1.arg1(), arg2), multiply(sum1.arg2(), arg2));
    Dout(dc::finish, result);
    return result;
  }
  else if (type2 == sumT)
  {
    Sum const& sum2 = static_cast<Sum const&>(arg2);
    Expression const& result = Sum::add(multiply(arg1, sum2.arg1()), multiply(arg1, sum2.arg2()));
    Dout(dc::finish, result);
    return result;
  }

  if (type1 == productT)
  {
    Product const* product1 = dynamic_cast<Product const*>(&arg1);
    if (type2 == productT)
    {
      Product const* product2 = dynamic_cast<Product const*>(&arg2);
      if (is_less(product1->arg1(), product2->arg1()))
      {
        Expression const& result = multiply(product1->arg1(), multiply(product1->arg2(), arg2));
        Dout(dc::finish, result);
        return result;
      }
      else if (is_less(product2->arg1(), product1->arg1()))
      {
        Expression const& result = multiply(product2->arg1(), multiply(arg1, product2->arg2()));
        Dout(dc::finish, result);
        return result;
      }
      else
      {
        Expression const& result = multiply(multiply(combine(product1->arg1(), product2->arg1()), product1->arg2()), product2->arg2());
        Dout(dc::finish, result);
        return result;
      }
    }
    if (is_less(product1->arg1(), arg2))
    {
      Expression const& result = multiply(product1->arg1(), multiply(product1->arg2(), arg2));
      Dout(dc::finish, result);
      return result;
    }
    else if (is_less(arg2, product1->arg1()))
    {
      Expression const& result = make_product(arg2, arg1);
      Dout(dc::finish, result);
      return result;
    }
    else
    {
      Expression const& result = multiply(combine(product1->arg1(), arg2), product1->arg2());
      Dout(dc::finish, result);
      return result;
    }
  }
  if (type2 == productT)
  {
    Product const* product2 = dynamic_cast<Product const*>(&arg2);
    if (is_less(arg1, product2->arg1()))
    {
      Expression const& result = make_product(arg1, arg2);
      Dout(dc::finish, result);
      return result;
    }
    else if (is_less(product2->arg1(), arg1))
    {
      Expression const& result = multiply(product2->arg1(), multiply(arg1, product2->arg2()));
      Dout(dc::finish, result);
      return result;
    }
    else
    {
      Expression const& result = multiply(combine(arg1, product2->arg1()), product2->arg2());
      Dout(dc::finish, result);
      return result;
    }
  }

  // Neither is a Sum or Product.
  if (is_less(arg1, arg2))
  {
    Expression const& result = make_product(arg1, arg2);
    Dout(dc::finish, result);
    return result;
  }
  else if (is_less(arg2, arg1))
  {
    Expression const& result = make_product(arg2, arg1);
    Dout(dc::finish, result);
    return result;
  }
  else if (Constant const* constant1 = dynamic_cast<Constant const*>(&arg1))
  {
    ASSERT(arg2.is_constant());
    Constant const& constant2 = static_cast<Constant const&>(arg2);
    Expression const& result = Constant::realize(constant1->enumerator_ * constant2.enumerator_, constant1->denominator_ * constant2.denominator_);
    Dout(dc::finish, result);
    return result;
  }
  else
  {
    Expression const& result = combine(arg1, arg2);
    Dout(dc::finish, result);
    return result;
  }
}

//static
Expression const& Product::negate(Expression const& arg)
{
  DoutEntering(dc::symbolic|continued_cf, "Product::negate(" << arg << ") --> ");

  Constant const& constant_factor = arg.get_constant_factor();
  if (Constant::is_one(constant_factor))
  {
    Expression const& result = multiply(Constant::s_cached_minus_one, arg);
    Dout(dc::finish, result);
    return result;
  }

  // If arg has a constant factor that is not one, then it must be a product (get_constant_factor() called on a Constant return one).
  ASSERT(arg.is_product());
  Product const& product = static_cast<Product const&>(arg);
  Constant const& negated_constant_factor = Constant::s_cached_minus_one * static_cast<Constant const&>(product.arg1());
  Expression const& result = make_product(negated_constant_factor, product.arg2());
  Dout(dc::finish, result);
  return result;
}

#ifdef SYMBOLIC_PRINTING
void Product::print_on(std::ostream& os) const
{
  if (Constant::is_minus_one(arg1_))
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
    bool need_parens = needs_parens(arg1_.precedence(), Precedence::product, before);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " * ";
    need_parens = needs_parens(arg2_.precedence(), Precedence::product, after);
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
}
#endif

} // namespace symbolic
