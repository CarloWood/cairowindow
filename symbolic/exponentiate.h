#pragma once

#include "expression_traits.h"
#include <utility>

namespace symbolic {

template<Expression Base, ConstantType Exponent>
struct exponentiate
{
  using type = std::remove_cvref_t<decltype(std::declval<Base>() ^ std::declval<Exponent>())>;
};

template<Expression Base, ConstantType Exponent>
using exponentiate_t = typename exponentiate<Base, Exponent>::type;

} // namespace symbolic
