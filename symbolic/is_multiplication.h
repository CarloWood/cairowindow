#pragma once

#include "Expression.h"
#include <type_traits>

namespace symbolic {

template<Expression E1, Expression E2>
class Multiplication;

template<typename T>
struct is_multiplication : std::false_type { };

template<typename T>
constexpr bool is_multiplication_v = is_multiplication<T>::value;

template<Expression E1, Expression E2>
struct is_multiplication<Multiplication<E1, E2>> : std::true_type { };

template<typename T>
concept MultiplicationType = is_multiplication_v<T>;

} // namespace symbolic
