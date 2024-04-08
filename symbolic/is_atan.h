#pragma once

#include "Expression.h"
#include <type_traits>

namespace symbolic {

template<Expression E>
class Atan;

template<typename T>
struct is_atan : std::false_type { };

template<typename T>
constexpr bool is_atan_v = is_atan<T>::value;

template<Expression E>
struct is_atan<Atan<E>> : std::true_type { };

template<typename T>
concept AtanType = is_atan_v<T>;

} // namespace symbolic
