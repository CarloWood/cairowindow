#pragma once

#include "is_power.h"

namespace symbolic {

template<typename T>
struct is_product : std::false_type { };

template<typename T>
constexpr bool is_product_v = is_product<T>::value;

template<Expression E1, Expression E2>
requires (((is_constant_v<E1> && !is_constant_zero_v<E1> && !is_constant_one_v<E1>) ||
           is_symbol_v<E1> || is_power_v<E1>) &&
          (is_symbol_v<E2> || is_power_v<E2> || is_product_v<E2>))
class Product;

template<Expression E1, Expression E2>
struct is_product<Product<E1, E2>> : std::true_type { };

template<typename E>
concept ProductType = is_product_v<E>;

template<typename T>
concept ProductLevelType = is_constant_v<T> || is_symbol_v<T> || is_power_v<T> || is_product_v<T>;

} // namespace symbolic
