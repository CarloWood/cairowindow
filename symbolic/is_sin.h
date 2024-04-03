#pragma once

#include "Expression.h"
#include <type_traits>

namespace symbolic {

template<Expression E>
class Sin;

template<typename T>
struct is_sin : std::false_type { };

template<typename T>
constexpr bool is_sin_v = is_sin<T>::value;

template<Expression E>
struct is_sin<Sin<E>> : std::true_type { };

template<typename T>
concept SinType = is_sin_v<T>;

} // namespace symbolic
