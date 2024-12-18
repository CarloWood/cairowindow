#pragma once

#include <string>

namespace enable_drawing {

class Function
{
 public:
  virtual double evaluate(double w) const = 0;
#ifdef CWDEBUG
  virtual std::string to_string() const = 0;
#endif

  // Provide call operator.
  double operator()(double w) const { return evaluate(w); }
};

} // namespace enable_drawing
