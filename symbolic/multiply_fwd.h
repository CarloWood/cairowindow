#pragma once

#include "expression_traits.h"
#include "Expression.h"

namespace symbolic {

enum IsProduct
{
  isProduct,            // Exclusively exists of (optionally one constant factor and) symbols as factors (Symbol and Power).
  not_a_Product         // All other multiplications.
};

template<Expression E1, Expression E2, IsProduct is_product =
  ((is_constant_v<E1> && !is_constant_v<E2>) || is_symbol_v<E1> || is_power_v<E1> || is_product_v<E1>) &&
  ((is_constant_v<E2> && !is_constant_v<E1>) || is_symbol_v<E2> || is_power_v<E2> || is_product_v<E2>) ? isProduct : not_a_Product>
struct multiply;

template<Expression E1, Expression E2>
using multiply_t = typename multiply<E1, E2>::type;

} // namespace symbolic
