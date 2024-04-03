#pragma once

#include <type_traits>

namespace symbolic {

template<int Id>
class Symbol;

template<typename T>
struct is_symbol : std::false_type { };

template<typename T>
constexpr bool is_symbol_v = is_symbol<T>::value;

template<int Id>
struct is_symbol<Symbol<Id>> : std::true_type { };

template<typename E>
concept SymbolType = is_symbol_v<E>;

} // namespace symbolic
