#pragma once

#include "UnaryOperator.h"
#include "debug.h"

namespace symbolic2 {

class Sin : public UnaryOperator<Sin>
{
 private:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Sin(Expression const& arg) : UnaryOperator<Sin>(arg) { }

  Precedence precedence() const override final { return Precedence::symbol; }

  ExpressionType type() const override final { return sinT; }

  double evaluate() const override final;
  Expression const& derivative(Symbol const& symbol) const override final;

#ifdef SYMBOLIC2_PRINTING
 public:
  void print_on(std::ostream& os) const override final
  {
    os << "sin(" << arg_ << ")";
  }
#endif
};

inline Expression const& sin(Expression const& arg)
{
  return Sin::realize(arg);
}

} // namespace symbolic2
