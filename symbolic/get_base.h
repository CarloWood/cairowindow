#pragma once

#include "Expression.h"
#include "is_constant.h"
#include "is_symbol.h"
#include "is_power.h"

namespace symbolic {

template<Expression E>
struct get_base
{
  using type = E;
};

template<Expression E>
using get_base_t = typename get_base<E>::type;

template<SymbolType Base, ConstantType Exponent>
struct get_base<Power<Base, Exponent>>
{
  using type = Base;
};

} // namespace symbolic
