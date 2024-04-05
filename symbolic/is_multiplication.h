#pragma once

#include "is_product.h"
#include "Expression.h"
#include <type_traits>

namespace symbolic {

template<typename T>
struct is_multiplication : std::false_type { };

template<typename T>
constexpr bool is_multiplication_v = is_multiplication<T>::value;

template<Expression E1, Expression E2>
requires ((!ProductLevelType<E1> || !ProductLevelType<E2>) && !is_multiplication_v<E1> && is_constant_factor_free_v<E2>)
class Multiplication;

template<int e, int d, Expression arg2>
struct is_constant_factor_free<Multiplication<Constant<e, d>, arg2>> : std::false_type { };

template<int e, int d, Expression E1, Expression arg2>
struct is_constant_factor_free<Multiplication<Product<Constant<e, d>, E1>, arg2>> : std::false_type { };

template<Expression E1, Expression E2>
struct is_multiplication<Multiplication<E1, E2>> : std::true_type { };

template<typename T>
concept MultiplicationType = is_multiplication_v<T>;

} // namespace symbolic
