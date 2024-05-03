#pragma once

#include "UnaryOperator.h"
#include "debug.h"

namespace symbolic {

class Exponential : public UnaryOperator<Exponential>
{
 private:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Exponential(Expression const& exponent) : UnaryOperator<Exponential>(exponent) { }

  Precedence precedence() const override final { return Precedence::power; }

  ExpressionType type() const override final { return expT; }

  double evaluate() const override final;
  Expression const& derivative(Symbol const& symbol) const override final;

#ifdef SYMBOLIC_PRINTING
 public:
  void print_on(std::ostream& os) const override final
  {
    os << "exp(" << arg_ << ")";
  }
#endif
};

inline Expression const& exp(Expression const& arg)
{
  return Exponential::realize(arg);
}

} // namespace symbolic
