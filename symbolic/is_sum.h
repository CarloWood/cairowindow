#pragma once

#include "is_constant.h"
#include "Expression.h"
#include <type_traits>

namespace symbolic {

template<typename T>
struct is_sum : std::false_type { };

template<typename T>
constexpr bool is_sum_v = is_sum<T>::value;

template<Expression E1, Expression E2>
requires (!is_sum_v<E1> && !is_constant_zero_v<E1> && !is_constant_v<E2>)
class Sum;

template<Expression E1, Expression E2>
struct is_sum<Sum<E1, E2>> : std::true_type { };

template<typename E>
concept SumType = is_sum_v<E>;

} // namespace symbolic
