#pragma once

#include <type_traits>

namespace symbolic {

template<typename T>
struct is_constant : std::false_type {};

template<typename T>
struct is_constant_zero : std::false_type {};

template<typename T>
struct is_constant_one : std::false_type {};

template<typename T>
struct is_constant_minus_one : std::false_type {};

template<typename T>
struct is_constant_less_than_zero : std::false_type {};

template<typename T>
constexpr bool is_constant_v = is_constant<T>::value;

template<typename T>
constexpr bool is_constant_zero_v = is_constant_zero<T>::value;

template<typename T>
constexpr bool is_constant_one_v = is_constant_one<T>::value;

template<typename T>
constexpr bool is_constant_minus_one_v = is_constant_minus_one<T>::value;

template<typename T>
constexpr bool is_constant_less_than_zero_v = is_constant_less_than_zero<T>::value;

// Forward declaration of Constant.
template<int Enumerator, int Denominator>
class Constant;

// Any Constant.
template<int Enumerator, int Denominator>
struct is_constant<Constant<Enumerator, Denominator>> : std::true_type {};

// Zero
template<>
struct is_constant_zero<Constant<0, 1>> : std::true_type {};

// One
template<>
struct is_constant_one<Constant<1, 1>> : std::true_type {};

// Minus one
template<>
struct is_constant_minus_one<Constant<-1, 1>> : std::true_type {};

// Less than zero
template<int Enumerator, int Denominator>
requires (Enumerator < 0)
struct is_constant_less_than_zero<Constant<Enumerator, Denominator>> : std::true_type {};

} // namespace symbolic
