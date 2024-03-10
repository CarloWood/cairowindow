#pragma once

namespace symbolic {

enum class precedence
{
  constant,
  symbol = constant,
  negation,
  product,
  division,
  exponentiation
};

enum before_or_after
{
  before,
  after
};

} // namespace symbolic
