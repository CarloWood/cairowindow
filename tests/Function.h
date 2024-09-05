#pragma once

#include <string>

namespace enable_drawing {

class Function
{
 public:
  virtual double operator()(double w) const = 0;
  virtual std::string to_string() const = 0;
};

} // namespace enable_drawing
