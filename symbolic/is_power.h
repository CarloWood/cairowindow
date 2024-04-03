#pragma once

#include "is_constant.h"
#include "is_symbol.h"
#include "Expression.h"

namespace symbolic {

template<SymbolType Base, ConstantType Exponent>
requires (!is_constant_zero_v<Exponent> && !is_constant_one_v<Exponent>)
class Power;

template<typename T>
struct is_power : std::false_type { };

template<typename T>
constexpr bool is_power_v = is_power<T>::value;

template<int Id, ConstantType Exponent>
struct is_power<Power<Symbol<Id>, Exponent>> : std::true_type { };

template<typename T>
concept PowerType = is_power_v<T>;

template<typename E>
concept SymbolPowerType = is_symbol_v<E> || is_power_v<E>;

} // namespace symbolic
