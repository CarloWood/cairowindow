#pragma once

#include "Product.h"
#include "debug.h"

namespace symbolic {

class Sum : public BinaryOperator<Sum>
{
 private:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Sum(Expression const& arg1, Expression const& arg2) : BinaryOperator<Sum>(arg1, arg2)
  {
    // First argument is not allowed to be zero.
    ASSERT(!Constant::is_zero(arg1));
    // Every Sum must have a non-Sum as first argument.
    ASSERT(!arg1.is_sum());
    // arg1 must always be "less than" arg2.
    ASSERT(Sum::is_less(arg1, arg2));
    // arg2 is not allowed to be a Sum that has a constant as first argument.
    ASSERT(!arg2.is_sum() || !arg2.arg1().is_constant());
  }

 private:
  Precedence precedence() const override final { return Precedence::sum; }

 public:
  static Expression const& add(Expression const& arg1, Expression const& arg2);
  static Expression const& combine(Expression const& arg1, Expression const& arg2);

  ExpressionType type() const override final { return sumT; }

  double evaluate() const override final
  {
    return arg1_.evaluate() + arg2_.evaluate();
  }

  Expression const& derivative(Symbol const& symbol) const override final;

#ifdef SYMBOLIC_PRINTING
 public:
  void print_on(std::ostream& os) const override final;
#endif
};

inline Expression const& operator+(Expression const& arg1, Expression const& arg2)
{
  return Sum::add(arg1, arg2);
}

inline Expression const& operator+(int arg1, Expression const& arg2)
{
  return Sum::add(Constant::realize(arg1), arg2);
}

inline Expression const& operator+(Expression const& arg1, int arg2)
{
  return Sum::add(arg1, Constant::realize(arg2));
}

inline Expression const& operator-(Expression const& arg1, Expression const& arg2)
{
  return Sum::add(arg1, Product::negate(arg2));
}

inline Expression const& operator-(int arg1, Expression const& arg2)
{
  return Sum::add(Constant::realize(arg1), Product::negate(arg2));
}

inline Expression const& operator-(Expression const& arg1, int arg2)
{
  return Sum::add(arg1, Constant::realize(-arg2));
}

} // namespace symbolic
