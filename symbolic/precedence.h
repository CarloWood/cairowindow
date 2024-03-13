#pragma once

namespace symbolic {

enum class precedence
{
  constant,
  symbol = constant,
  negation,
  exponentiation,
  product,
  division
};

enum before_or_after
{
  before,
  after
};

} // namespace symbolic
