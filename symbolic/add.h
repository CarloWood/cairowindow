#pragma once

#include "Constant.h"
#include "is_less_Sum.h"
#include "add_equals.h"
#include "Sum.h"

namespace symbolic {

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

} // namespace symbolic
