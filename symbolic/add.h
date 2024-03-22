#pragma once

#include "expression_traits.h"
#include "Constant.h"

namespace symbolic {

template<Expression E1, Expression E2>
struct add_equals;

template<Expression E1, Expression E2>
using add_equals_t = typename add_equals<E1, E2>::type;

// This default is for when E1 nor E2 is a Sum. If one or both are a Sum then
// one of the specializations below will be used.
template<Expression E1, Expression E2>
class add
{
  static consteval auto eval()
  {
    if constexpr (is_less_Sum_v<E1, E2>)
    {
      if constexpr (is_constant_zero_v<E1>)
        return E2::instance();
      else
        return Sum<E1, E2>::instance();
    }
    else if constexpr (is_less_Sum_v<E2, E1>)
    {
      if constexpr (is_constant_zero_v<E2>)
        return E1::instance();
      else
        return Sum<E2, E1>::instance();
    }
    else
      return add_equals_t<E1, E2>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2>
using add_t = typename add<E1, E2>::type;

template<Expression E1, Expression E2, Expression E3>
class add<E1, Sum<E2, E3>>
{
  static consteval auto eval()
  {
    if constexpr (is_less_Sum_v<E1, E2>)
    {
      if constexpr (is_constant_zero_v<E1>)
        return Sum<E2, E3>{};
      else
        return Sum<E1, Sum<E2, E3>>{};
    }
    else if constexpr (is_less_Sum_v<E2, E1>)
      return Sum<E2, add_t<E1, E3>>{};
    else
      return add_t<add_equals_t<E1, E2>, E3>{};
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2, Expression E3>
class add<Sum<E1, E2>, E3>
{
  static consteval auto eval()
  {
    if constexpr (is_less_Sum_v<E1, E3>)
      return Sum<E1, add_t<E3, E2>>::instance();
    else if constexpr (is_less_Sum_v<E3, E1>)
      return Sum<E3, Sum<E1, E2>>::instance();
    else
      return add_t<add_equals_t<E1, E3>, E2>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2, Expression E3, Expression E4>
class add<Sum<E1, E2>, Sum<E3, E4>>
{
  static consteval auto eval()
  {
    if constexpr (is_less_Sum_v<E1, E3>)
      return Sum<E1, add_t<E2, Sum<E3, E4>>>{};
    else if constexpr (is_less_Sum_v<E3, E1>)
      return Sum<E3, add_t<Sum<E1, E2>, E4>>{};
    else
      return add_t<add_t<add_equals_t<E1, E3>, E2>, E4>{};
  }

 public:
  using type = decltype(eval());
};

template<Expression E>
struct get_constant_factor
{
  using type = Constant<1, 1>;
};

template<Expression E>
using get_constant_factor_t = typename get_constant_factor<E>::type;

template<int Enumerator1, int Denominator1, Expression E2>
struct get_constant_factor<Product<Constant<Enumerator1, Denominator1>, E2>>
{
  using type = Constant<Enumerator1, Denominator1>;
};

template<Expression E>
struct get_nonconstant_factor
{
  using type = E;
};

template<Expression E>
using get_nonconstant_factor_t = typename get_nonconstant_factor<E>::type;

template<int Enumerator1, int Denominator1, Expression E2>
struct get_nonconstant_factor<Product<Constant<Enumerator1, Denominator1>, E2>>
{
  using type = E2;
};

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
    else
      return Product<ConstantFactor, non_constant_factor1>{};
  }

 public:
  using type = decltype(eval());
};

} // namespace symbolic
