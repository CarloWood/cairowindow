#pragma once

#include "BinaryOperator.h"

namespace symbolic {

class Product : public BinaryOperator<Product>
{
 private:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Product(Expression const& arg1, Expression const& arg2) : BinaryOperator<Product>(arg1, arg2)
  {
    // First argument is not allowed to be zero or one.
    ASSERT(!Constant::is_zero(arg1) && !Constant::is_one(arg1));
    // Every Product must have a non-Product as first argument.
    ASSERT(!arg1.is_product());
    // arg1 must always be "less than" arg2.
    ASSERT(Product::is_less(arg1, arg2));
    // arg2 is not allowed to be a Product that has a constant as first argument.
    ASSERT(!arg2.is_product() || !arg2.arg1().is_constant());
  }

 private:
  Precedence precedence() const override final
  {
    Constant const* constant_factor = dynamic_cast<Constant const*>(&arg1_);
    return (constant_factor && constant_factor->is_negative()) ? Precedence::negation : Precedence::product;
  }

 public:
  static Expression const& combine(Expression const& arg1, Expression const& arg2);
  static Expression const& multiply(Expression const& arg1, Expression const& arg2);
  static Expression const& negate(Expression const& arg);

  static Expression const& make_product(Expression const& arg1, Expression const& arg2);

  ExpressionType type() const override final { return productT; }

  Expression const& get_nonconstant_factor() const override final
  {
    return arg1_.is_constant() ? arg2_ : *this;
  }

  Constant const& get_constant_factor() const override final
  {
    Constant const* constant_factor = dynamic_cast<Constant const*>(&arg1_);
    return constant_factor ? *constant_factor : Constant::s_cached_one;
  }

  static bool has_constant_factor(Expression const& arg)
  {
    Product const* product = dynamic_cast<Product const*>(&arg);
    return product && product->arg1().is_constant();
  }

  double evaluate() const override final
  {
    return arg1_.evaluate() * arg2_.evaluate();
  }

  Expression const& derivative(Symbol const& symbol) const override final;

#ifdef SYMBOLIC_PRINTING
 public:
  void print_on(std::ostream& os) const override final;
#endif
};

inline Expression const& operator*(Expression const& arg1, Expression const& arg2)
{
  return Product::multiply(arg1, arg2);
}

inline Expression const& operator*(int arg1, Expression const& arg2)
{
  return Product::multiply(Constant::realize(arg1), arg2);
}

inline Expression const& operator*(Expression const& arg1, int arg2)
{
  return Product::multiply(arg1, Constant::realize(arg2));
}

inline Expression const& operator/(Expression const& arg1, Expression const& arg2)
{
  return Product::multiply(arg1, arg2^-1);
}

inline Expression const& operator/(int arg1, Expression const& arg2)
{
  return Product::multiply(Constant::realize(arg1), arg2^-1);
}

inline Expression const& operator/(Expression const& arg1, int arg2)
{
  return Product::multiply(arg1, Constant::realize(1, arg2));
}

inline Expression const& operator-(Expression const& arg)
{
  return Product::negate(arg);
}

} // namespace symbolic
