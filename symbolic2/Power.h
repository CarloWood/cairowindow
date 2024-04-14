#pragma once

#include "BinaryOperator.h"
#include <cmath>

namespace symbolic2 {

class Power : public BinaryOperator<Power>
{
 private:
  Precedence precedence() const override final { return Precedence::power; }

  // When it is known that exponent isn't zero or one.
  static Expression const& make_power2(Expression const& base, Constant const& exponent);

 public:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Power(Expression const& arg1, Expression const& arg2);

  static Expression const& make_power(Expression const& base, Constant const& exponent);

  ExpressionType type() const override final { return powerT; }

  Expression const& get_base() const override final { return arg1_; }
  Constant const& get_exponent() const override final { return static_cast<Constant const&>(arg2_); }

  double evaluate() const override final
  {
    return std::pow(arg1_.evaluate(), arg2_.evaluate());
  }

  Expression const& differentiate(Symbol const& symbol) const override final;

#ifdef SYMBOLIC2_PRINTING
 public:
  void print_on(std::ostream& os) const override final;
#endif
};

inline Expression const& operator^(Expression const& arg1, Expression const& arg2)
{
  Constant const* base = dynamic_cast<Constant const*>(&arg1);
  Constant const* exponent = dynamic_cast<Constant const*>(&arg2);

  // Can only raise to a constant power.
  ASSERT(exponent);

  if (!base)
    return Power::make_power(arg1, *exponent);

  // For example, when trying to calculate: (27/8)^(5/3)

  int e1 = base->enumerator_;           // 27
  int d1 = base->denominator_;          // 8
  int e2 = exponent->enumerator_;       // 5
  int d2 = exponent->denominator_;      // 3

  if (e2 < 0)
  {
    std::swap(e1, d1);
    e2 = -e2;
  }

  long e_nb = std::lround(std::pow(static_cast<double>(e1), 1.0 / d2));         // 27^(1/3) = 3
  long d_nb = std::lround(std::pow(static_cast<double>(d1), 1.0 / d2));         // 8^(1/3) = 2

  if (std::pow(e_nb, d2) != e1 ||       // 3^3 == 27
      std::pow(d_nb, d2) != d1)         // 2^3 == 8
    DoutFatal(dc::core, "Can't exponentiate " << arg1 << " to the power " << arg2);

  return Constant::realize(std::pow(e_nb, e2), std::pow(d_nb, e2));     // (3^5) / (2^5) = 243/32
}

inline Expression const& operator^(Expression const& arg1, int exponent)
{
  return arg1 ^ Constant::realize(exponent);
}

} // namespace symbolic2

namespace utils {

inline symbolic2::Expression const& square(symbolic2::Expression const& expression)
{
  static symbolic2::Constant const& cached_two = symbolic2::Constant::realize(2);
  return symbolic2::Power::make_power(expression, cached_two);
}

}
