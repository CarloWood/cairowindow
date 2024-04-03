#pragma once

#include "Expression.h"
#include "Constant.h"

namespace symbolic {

template<Expression E>
struct get_exponent
{
  using type = Constant<1, 1>;
};

template<Expression E>
using get_exponent_t = typename get_exponent<E>::type;

} // namespace symbolic
