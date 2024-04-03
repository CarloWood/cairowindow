#pragma once

#include "Expression.h"
#include <type_traits>

namespace symbolic {

template<Expression E>
class Cos;

template<typename T>
struct is_cos : std::false_type { };

template<typename T>
constexpr bool is_cos_v = is_cos<T>::value;

template<Expression E>
struct is_cos<Cos<E>> : std::true_type { };

template<typename T>
concept CosType = is_cos_v<T>;

} // namespace symbolic
