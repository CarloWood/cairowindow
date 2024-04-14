#pragma once

#include "UnaryOperator.h"
#include "debug.h"

namespace symbolic2 {

class Log : public UnaryOperator<Log>
{
 private:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Log(Expression const& arg) : UnaryOperator<Log>(arg) { }

  Precedence precedence() const override final { return Precedence::symbol; }

  ExpressionType type() const override final { return logT; }

  double evaluate() const override final;
  Expression const& differentiate(Symbol const& symbol) const override final;

#ifdef SYMBOLIC2_PRINTING
 public:
  void print_on(std::ostream& os) const override final
  {
    os << "log(" << arg_ << ")";
  }
#endif
};

inline Expression const& log(Expression const& arg)
{
  return Log::realize(arg);
}

} // namespace symbolic2
