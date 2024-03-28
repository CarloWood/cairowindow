#pragma once

#include "is_constant.h"
#include "Expression.h"

namespace symbolic {

template<typename E>
concept ConstantType = is_constant_v<std::remove_const_t<E>>;

template<int Id>
class Symbol;

template<Expression Base, ConstantType Exponent>
requires (!is_constant_zero_v<Exponent> && !is_constant_one_v<Exponent>)
class Power;

template<Expression E1, Expression E2>
class Sum;

template<Expression Base, ConstantType Exponent>
requires (!is_constant_zero_v<Exponent> && !is_constant_one_v<Exponent>)
class Exponentiation;

template<typename T>
struct is_symbol : std::false_type { };

template<typename T>
constexpr bool is_symbol_v = is_symbol<T>::value;

template<int Id>
struct is_symbol<Symbol<Id>> : std::true_type { };

template<typename E>
concept SymbolType = is_symbol_v<E>;

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

template<typename T>
struct is_sum : std::false_type { };

template<typename T>
constexpr bool is_sum_v = is_sum<T>::value;

template<Expression E1, Expression E2>
struct is_sum<Sum<E1, E2>> : std::true_type { };

template<typename E>
concept SumType = is_sum_v<E>;

template<typename T>
struct is_exponentiation : std::false_type { };

template<typename T>
constexpr bool is_exponentiation_v = is_exponentiation<T>::value;

template<Expression Base, ConstantType Exponent>
struct is_exponentiation<Exponentiation<Base, Exponent>> : std::true_type { };

template<typename T>
concept ExponentiationType = is_exponentiation_v<T>;

template<Expression E1, Expression E2>
struct is_less_Sum : std::false_type { };

template<Expression E1, Expression E2>
constexpr bool is_less_Sum_v = is_less_Sum<E1, E2>::value;

// A Constant is (only) less than a non-constant.
template<int Enumerator1, int Denominator1, Expression E2>
requires (!is_constant_v<E2>)
struct is_less_Sum<Constant<Enumerator1, Denominator1>, E2> : std::true_type { };

// Compare two symbols.
template<int Id1, int Id2>
requires (Id1 < Id2)
struct is_less_Sum<Symbol<Id1>, Symbol<Id2>> : std::true_type { };

// Compare a symbol with something not a constant or symbol.
template<int Id1, Expression E2>
requires (!is_constant_v<E2> && !is_symbol_v<E2> && !is_product_v<E2>)
struct is_less_Sum<Symbol<Id1>, E2> : std::true_type { };

template<int Id1, Expression E2, Expression E3>
requires (!is_constant_v<E2> || is_less_Sum_v<Symbol<Id1>, E3>)
struct is_less_Sum<Symbol<Id1>, Product<E2, E3>> : std::true_type { };

template<Expression E1, Expression E2, int Id3>
requires (is_constant_v<E1> && is_less_Sum_v<E2, Symbol<Id3>>)
struct is_less_Sum<Product<E1, E2>, Symbol<Id3>> : std::true_type { };

template<int Id1, ConstantType Exponent1, int Id2, ConstantType Exponent2>
requires (is_less_Sum_v<Symbol<Id1>, Symbol<Id2>> || (!is_less_Sum_v<Symbol<Id2>, Symbol<Id1>> &&
      Exponent1::s_enumerator * Exponent2::s_denominator < Exponent2::s_enumerator * Exponent1::s_denominator))
struct is_less_Sum<Power<Symbol<Id1>, Exponent1>, Power<Symbol<Id2>, Exponent2>> : std::true_type { };

template<int Id1, ConstantType Exponent1, Expression E2>
requires (!is_constant_v<E2> && !is_symbol_v<E2> && !is_power_v<E2> && !is_product_v<E2>)
struct is_less_Sum<Power<Symbol<Id1>, Exponent1>, E2> : std::true_type { };

template<int Id1, ConstantType Exponent1, Expression E2, Expression E3>
requires (!is_constant_v<E2> || is_less_Sum_v<Power<Symbol<Id1>, Exponent1>, E3>)
struct is_less_Sum<Power<Symbol<Id1>, Exponent1>, Product<E2, E3>> : std::true_type { };

template<Expression E1, Expression E2, int Id, ConstantType Exponent>
requires (is_constant_v<E1> && is_less_Sum_v<E2, Power<Symbol<Id>, Exponent>>)
struct is_less_Sum<Product<E1, E2>, Power<Symbol<Id>, Exponent>> : std::true_type { };

template<Expression E1, Expression E2, Expression E3, Expression E4>
requires (
    (!is_constant_v<E1> && !is_constant_v<E3> &&
     (is_less_Sum_v<E1, E3> || (!is_less_Sum_v<E3, E1> && is_less_Sum_v<E2, E4>))) ||
    (!is_constant_v<E1> && is_constant_v<E3> && is_less_Sum_v<Product<E1, E2>, E4>) ||
    (is_constant_v<E1> && !is_constant_v<E3> && is_less_Sum_v<E2, Product<E3, E4>>) ||
    (is_constant_v<E1> && is_constant_v<E3> && is_less_Sum_v<E2, E4>))
struct is_less_Sum<Product<E1, E2>, Product<E3, E4>> : std::true_type { };

template<Expression E1, Expression E2, Expression E3>
requires (!is_constant_v<E3> && !is_symbol_v<E3> && !is_power_v<E3> && !is_product_v<E3>)
struct is_less_Sum<Product<E1, E2>, E3> : std::true_type { };

template<Expression E1, Expression E2, Expression E3, Expression E4>
requires (is_less_Sum_v<E1, E3> || (!is_less_Sum_v<E3, E1> && is_less_Sum_v<E2, E4>))
struct is_less_Sum<Sum<E1, E2>, Sum<E3, E4>> : std::true_type { };

template<Expression E1, Expression E2, Expression E3>
requires (!is_constant_v<E3> && !is_symbol_v<E3> && !is_power_v<E3> && !is_product_v<E3>)
struct is_less_Sum<Sum<E1, E2>, E3> : std::true_type { };

template<Expression E1, Expression E2>
struct is_less_Product : std::false_type { };

template<Expression E1, Expression E2>
constexpr bool is_less_Product_v = is_less_Product<E1, E2>::value;

template<int Enumerator1, int Denominator1, Expression E2>
requires (!is_constant_v<E2>)
struct is_less_Product<Constant<Enumerator1, Denominator1>, E2> : std::true_type { };

template<int Id1, int Id2>
requires (Id1 < Id2)
struct is_less_Product<Symbol<Id1>, Symbol<Id2>> : std::true_type { };

template<int Id1, int Id2, ConstantType Exponent>
requires (Id1 < Id2)
struct is_less_Product<Symbol<Id1>, Power<Symbol<Id2>, Exponent>> : std::true_type { };

template<int Id1, ConstantType Exponent1, int Id2>
requires (Id1 < Id2)
struct is_less_Product<Power<Symbol<Id1>, Exponent1>, Symbol<Id2>> : std::true_type { };

template<int Id1, ConstantType Exponent1, int Id2, ConstantType Exponent2>
requires (Id1 < Id2)
struct is_less_Product<Power<Symbol<Id1>, Exponent1>, Power<Symbol<Id2>, Exponent2>> : std::true_type { };

template<Expression E>
struct get_exponent;

template<Expression E>
using get_exponent_t = typename get_exponent<E>::type;

template<SymbolPowerType E>
struct get_base;

template<SymbolPowerType E>
using get_base_t = typename get_base<E>::type;

} // namespace symbolic
