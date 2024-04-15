#include "sys.h"
#include "Product.h"
#include "Sum.h"
#include "Constant.h"
#include "Symbol.h"
#include "debug.h"

namespace symbolic2 {

//static
Expression const& Sum::combine(Expression const& arg1, Expression const& arg2)
{
  DoutEntering(dc::symbolic|continued_cf, "Sum::combine(" << arg1 << ", " << arg2 << ") --> ");

  Constant const* constant1 = dynamic_cast<Constant const*>(&arg1);
  if (constant1)
  {
    ASSERT(arg2.type() == constantT);
    Expression const& result = *constant1 + static_cast<Constant const&>(arg2);
    Dout(dc::finish, result);
    return result;
  }

  Expression const& factor1 = arg1.get_nonconstant_factor();
  Expression const& factor2 = arg2.get_nonconstant_factor();

  // Only call Sum::combine for two expression that have the same non-constant expression as factor.
  ASSERT(&factor1 == &factor2);

  Expression const& result = Product::multiply(arg1.get_constant_factor() + arg2.get_constant_factor(), factor1);
  Dout(dc::finish, result);
  return result;
}

//static
Expression const& Sum::add(Expression const& arg1, Expression const& arg2)
{
  DoutEntering(dc::symbolic|continued_cf, "Sum::add(" << arg1 << ", " << arg2 << ") --> ");

  Sum const* sum1 = dynamic_cast<Sum const*>(&arg1);
  Sum const* sum2 = dynamic_cast<Sum const*>(&arg2);

  if (!sum1)
  {
    if (!sum2)
    {
      // If neither argument is a Sum, then we have to create one
      // where the first argument is less than the second.
      if (is_less(arg1, arg2))
      {
        if (Constant::is_zero(arg1))
        {
          Dout(dc::finish, arg2);
          return arg2;
        }
        else
        {
          Expression const& result = realize(arg1, arg2);
          Dout(dc::finish, result);
          return result;
        }
      }
      else if (is_less(arg2, arg1))
      {
        if (Constant::is_zero(arg2))
        {
          Dout(dc::finish, arg1);
          return arg1;
        }
        else
        {
          Expression const& result = realize(arg2, arg1);
          Dout(dc::finish, result);
          return result;
        }
      }
      else
      {
        // If arg1 and arg2 are equal (up to a constant factor) they can be combined:
        Expression const& result = combine(arg1, arg2);
        Dout(dc::finish, result);
        return result;
      }
    }
    else
    {
      if (is_less(arg1, sum2->arg1_))
      {
        if (Constant::is_zero(arg1))
        {
          Dout(dc::finish, arg2);
          return arg2;
        }
        else
        {
          Expression const& result = realize(arg1, *sum2);
          Dout(dc::finish, result);
          return result;
        }
      }
      else if (is_less(sum2->arg1_, arg1))
      {
        Expression const& result = realize(sum2->arg1_, add(arg1, sum2->arg2_));
        Dout(dc::finish, result);
        return result;
      }
      else
      {
        Expression const& result = add(combine(arg1, sum2->arg1_), sum2->arg2_);
        Dout(dc::finish, result);
        return result;
      }
    }
  }
  else if (sum2)
  {
    if (is_less(sum1->arg1_, sum2->arg1_))
    {
      Expression const& result = add(sum1->arg1_, add(sum1->arg2_, *sum2));
      Dout(dc::finish, result);
      return result;
    }
    else if (is_less(sum2->arg1_, sum1->arg1_))
    {
      Expression const& result = add(sum2->arg1_, add(*sum1, sum2->arg2_));
      Dout(dc::finish, result);
      return result;
    }
    else
    {
      Expression const& result = add(add(combine(sum1->arg1_, sum2->arg1_), sum1->arg2_), sum2->arg2_);
      Dout(dc::finish, result);
      return result;
    }
  }
  else
  {
    if (is_less(sum1->arg1_, arg2))
    {
      Expression const& result = realize(sum1->arg1_, add(arg2, sum1->arg2_));
      Dout(dc::finish, result);
      return result;
    }
    else if (is_less(arg2, sum1->arg1_))
    {
      Expression const& result = realize(arg2, *sum1);
      Dout(dc::finish, result);
      return result;
    }
    else
    {
      Expression const& result = add(combine(sum1->arg1_, arg2), sum1->arg2_);
      Dout(dc::finish, result);
      return result;
    }
  }
}

Expression const& Sum::derivative(Symbol const& symbol) const
{
  return add(arg1_.derivative(symbol), arg2_.derivative(symbol));
}

#ifdef SYMBOLIC2_PRINTING
void Sum::print_on(std::ostream& os) const
{
  bool need_parens = needs_parens(arg1_.precedence(), Precedence::sum, before);
  if (need_parens)
    os << '(';
  arg1_.print_on(os);
  if (need_parens)
    os << ')';
  Constant const& constant_factor = arg2_.get_constant_factor();
  if (constant_factor.is_negative())
  {
    os << " - ";
    Expression const& term = Product::negate(arg2_);
    need_parens = needs_parens(term.precedence(), Precedence::difference, after);
    if (need_parens)
      os << '(';
    term.print_on(os);
    if (need_parens)
      os << ')';
    return;
  }
  else if (arg2_.is_sum())
  {
    Constant const& constant_factor = arg2_.arg1().get_constant_factor();
    if (constant_factor.is_negative())
    {
      os << " - ";
      Expression const& sum_with_negated_first_term = Sum::add(Product::negate(arg2_.arg1()), arg2_.arg2());
      need_parens = needs_parens(sum_with_negated_first_term.precedence(), Precedence::sum, after);
      if (need_parens)
        os << '(';
      sum_with_negated_first_term.print_on(os);
      if (need_parens)
        os << ')';
      return;
    }
  }
  os << " + ";
  need_parens = needs_parens(arg2_.precedence(), Precedence::sum, after);
  if (need_parens)
    os << '(';
  arg2_.print_on(os);
  if (need_parens)
    os << ')';
}
#endif

} // namespace symbolic2
