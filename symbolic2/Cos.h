#pragma once

#include "UnaryOperator.h"
#include "debug.h"

namespace symbolic2 {

class Cos : public UnaryOperator<Cos>
{
 private:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Cos(Expression const& arg) : UnaryOperator<Cos>(arg) { }

  Precedence precedence() const override final { return Precedence::symbol; }

  ExpressionType type() const override final { return cosT; }

  double evaluate() const override final;
  Expression const& differentiate(Symbol const& symbol) const override final;

#ifdef SYMBOLIC2_PRINTING
 public:
  void print_on(std::ostream& os) const override final
  {
    os << "cos(" << arg_ << ")";
  }
#endif
};

inline Expression const& cos(Expression const& arg)
{
  return Cos::realize(arg);
}

} // namespace symbolic2
