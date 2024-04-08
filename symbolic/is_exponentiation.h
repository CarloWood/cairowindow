#pragma once

#include "is_multiplication.h"
#include "Expression.h"

namespace symbolic {

template<Expression Base, ConstantType Exponent>
requires (!ProductLevelType<Base> && !is_constant_zero_v<Exponent> && !is_constant_one_v<Exponent> && !is_multiplication_v<Base>)
class Exponentiation;

template<typename T>
struct is_exponentiation : std::false_type { };

template<typename T>
constexpr bool is_exponentiation_v = is_exponentiation<T>::value;

template<Expression Base, ConstantType Exponent>
struct is_exponentiation<Exponentiation<Base, Exponent>> : std::true_type { };

template<typename T>
concept ExponentiationType = is_exponentiation_v<T>;

} // namespace symbolic
