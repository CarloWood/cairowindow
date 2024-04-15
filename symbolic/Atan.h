#pragma once

#include "UnaryOperator.h"
#include "debug.h"

namespace symbolic {

class Atan : public UnaryOperator<Atan>
{
 private:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Atan(Expression const& arg) : UnaryOperator<Atan>(arg) { }

  Precedence precedence() const override final { return Precedence::symbol; }

  ExpressionType type() const override final { return atanT; }

  double evaluate() const override final;
  Expression const& derivative(Symbol const& symbol) const override final;

#ifdef SYMBOLIC_PRINTING
 public:
  void print_on(std::ostream& os) const override final
  {
    os << "atan(" << arg_ << ")";
  }
#endif
};

inline Expression const& atan(Expression const& arg)
{
  return Atan::realize(arg);
}

} // namespace symbolic
