#pragma once

#include "is_constant.h"
#include "is_sum.h"
#include "get_constant_factor.h"
#include "get_nonconstant_factor.h"

namespace symbolic {

template<Expression E1, Expression E2>
class add_equals
{
  using constant_factor1 = get_constant_factor_t<E1>;
  using non_constant_factor1 = get_nonconstant_factor_t<E1>;
  using constant_factor2 = get_constant_factor_t<E2>;
  using non_constant_factor2 = get_nonconstant_factor_t<E2>;
  using ConstantFactor = typename add<constant_factor1, constant_factor2>::type;

  static_assert(!is_sum_v<E1>, "The first term of add_equals must not be a Sum.");
  static_assert(!is_constant_v<E1> || is_constant_v<E2>, "Equals means that if E1 is a constant then so should be E2!");
  static_assert(is_constant_v<E1> || std::is_same_v<non_constant_factor1, non_constant_factor2>, "Expected the same non-constant type.");

  static consteval auto eval()
  {
    if constexpr (is_constant_v<E1>)
      return add<E1, E2>::type::instance();
    else if constexpr (is_constant_zero_v<ConstantFactor>)
      return Constant<0, 1>::instance();
    else if constexpr (is_constant_one_v<ConstantFactor>)
      return non_constant_factor1::instance();
    else if constexpr (is_symbol_v<non_constant_factor1> || is_power_v<non_constant_factor1> || is_product_v<non_constant_factor1>)
      return Product<ConstantFactor, non_constant_factor1>::instance();
    else
      return Multiplication<ConstantFactor, non_constant_factor1>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2>
using add_equals_t = typename add_equals<E1, E2>::type;

} // namespace symbolic
