#pragma once

#include "is_product.h"
#include "is_exponentiation.h"
#include "is_multiplication.h"
#include "is_sum.h"
#include "is_sin.h"
#include "is_cos.h"
#include "Constant.h"
#include "Symbol.h"
#include "Power.h"
#include "Product.h"
#include "Sum.h"
#include "Multiplication.h"
#include "Exponentiation.h"
#include "functions.h"

namespace symbolic {

// Two expressions are only "equal" when they differ at most a constant factor.
// Any other ordering is therefore arbitrary - as long as terms that only
// differ a constant factor are brought together.
//
// Neverless, we use the following ordering:
//
//   Constant < anything
//   Symbol < Power < get_nonconstant_factor_t<Product> < Exponentiation < get_nonconstant_factor_t<Multiplication> < Sin < Cos < Sum.
//
// where we can assume that a Power and Exponentiation never have an exponent of 0 or 1,
// and the Product nor Multiplication begin with a constant.
//
// Note that Product < Multiplication in all cases, because it is impossible
// for the Multiplication to exist of a Constant times something that is less
// than a Product.
//
// In case two arguments are the same, the following holds:
//
//   1) two Constant's : always equal (aka, compare by "Id" of -1).
//   2) two Symbol's : only equal when their Id is the same (aka, compare as Id1 < Id2).
//   3) two Power's : only equal when exactly equal: first compare the Base, then the exponent.
//   4) two Product's : equal when the non-constant factor is equal.
//   5) two Exponentiation's : only equal when exactly equal: first compare the Base, then the exponent.
//   6) two Multiplication's : equal when the non-constant factor is equal.
//   7) two Sin : only equal when exactly equal (the arguments must exactly equal).
//   8) two Cos : only equal when exactly equal (the arguments must exactly equal).
//
template<int v>
struct ValueType
{
  static constexpr int value = v;
};

template<Expression E1, Expression E2, bool constants_compare_equal>
struct is_less_Sum : /*std::false_type*/ ValueType<0> { };

template<Expression E1, Expression E2>
constexpr bool is_less_Sum_v = is_less_Sum<E1, E2, true>::value;

// A Constant is (only) less than a non-constant.
template<int Enumerator1, int Denominator1, Expression E2, bool constants_compare_equal>
requires (!is_constant_v<E2>)
struct is_less_Sum<Constant<Enumerator1, Denominator1>, E2, constants_compare_equal> : /*std::true_type*/ ValueType<1> { };

// Unless Constant's need to be compared too (exact comparison).
template<int Enumerator1, int Denominator1, int Enumerator2, int Denominator2>
requires (Enumerator1 * Denominator2 < Enumerator2 * Denominator1)
struct is_less_Sum<Constant<Enumerator1, Denominator1>, Constant<Enumerator2, Denominator2>, false> : /*std::true_type*/ ValueType<2> { };

// A Symbol is always less than anything that isn't a Constant, Symbol, Product or Multiplication.
template<int Id1, Expression E2, bool constants_compare_equal>
requires (!is_constant_v<E2> && !is_symbol_v<E2> && !is_product_v<E2> && !is_multiplication_v<E2>)
struct is_less_Sum<Symbol<Id1>, E2, constants_compare_equal> : /*std::true_type*/ ValueType<3> { };

// A Symbol is only less than another Symbol if its Id is less.
template<int Id1, int Id2, bool constants_compare_equal>
requires (Id1 < Id2)
struct is_less_Sum<Symbol<Id1>, Symbol<Id2>, constants_compare_equal> : /*std::true_type*/ ValueType<4> { };

// Compare a Symbol with a Product: always less if the Product doesn't have a constant factor;
// otherwise when it is less than the non-constant factor.
template<int Id1, Expression E2, Expression E3, bool constants_compare_equal>
requires (!is_constant_v<E2> || is_less_Sum_v<Symbol<Id1>, E3>)
struct is_less_Sum<Symbol<Id1>, Product<E2, E3>, constants_compare_equal> : /*std::true_type*/ ValueType<5> { };

// Compare a Symbol with a Multiplication: always less if the Multiplication doesn't have a constant factor;
// otherwise when it is less than the non-constant factor.
template<int Id1, Expression E2, Expression E3, bool constants_compare_equal>
requires (!is_constant_v<E2>)
struct is_less_Sum<Symbol<Id1>, Multiplication<E2, E3>, constants_compare_equal> : /*std::true_type*/ ValueType<6> { };

// A Power is always less than anything that isn't a Constant, Symbol, Power, Product or Multiplication.
template<int Id1, ConstantType Exponent, Expression E2, bool constants_compare_equal>
requires (!is_constant_v<E2> && !is_symbol_v<E2> && !is_power_v<E2> && !is_product_v<E2> && !is_multiplication_v<E2>)
struct is_less_Sum<Power<Symbol<Id1>, Exponent>, E2, constants_compare_equal> : /*std::true_type*/ ValueType<7> { };

// A Power is only less than another Power if its Id is less, or
// its Id is equal and its exponent is less.
template<int Id1, ConstantType Exponent1, int Id2, ConstantType Exponent2, bool constants_compare_equal>
requires (Id1 < Id2 || (Id1 == Id2 && is_less_Sum<Exponent1, Exponent2, false>::value != 0))
struct is_less_Sum<Power<Symbol<Id1>, Exponent1>, Power<Symbol<Id2>, Exponent2>, constants_compare_equal> : /*std::true_type*/ ValueType<8> { };

// Compare a Power with a Product: always less if the Product doesn't have a constant factor;
// otherwise when it is less than the non-constant factor.
template<int Id1, ConstantType Exponent, Expression E2, Expression E3, bool constants_compare_equal>
requires (!is_constant_v<E2> || is_less_Sum_v<Power<Symbol<Id1>, Exponent>, E3>)
struct is_less_Sum<Power<Symbol<Id1>, Exponent>, Product<E2, E3>, constants_compare_equal> : /*std::true_type*/ ValueType<9> { };

// Compare a Power with a Multiplication: always less if the Multiplication doesn't have a constant factor;
// otherwise when it is less than the non-constant factor.
template<int Id1, ConstantType Exponent, Expression E2, Expression E3, bool constants_compare_equal>
requires (!is_constant_v<E2>)
struct is_less_Sum<Power<Symbol<Id1>, Exponent>, Multiplication<E2, E3>, constants_compare_equal> : /*std::true_type*/ ValueType<10> { };

// A Product is always less than anything that isn't a Constant, Symbol, Power or Product.
template<Expression E1, Expression E2, Expression E3, bool constants_compare_equal>
requires (!is_constant_v<E3> && !is_symbol_v<E3> && !is_power_v<E3> && !is_product_v<E3>)
struct is_less_Sum<Product<E1, E2>, E3, constants_compare_equal> : /*std::true_type*/ ValueType<11> { };

// Compare Product with Symbol.
template<Expression E1, Expression E2, int Id3, bool constants_compare_equal>
requires (is_constant_v<E1> && is_less_Sum_v<E2, Symbol<Id3>>)
struct is_less_Sum<Product<E1, E2>, Symbol<Id3>, constants_compare_equal> : /*std::true_type*/ ValueType<12> { };

// Compare Product with Power.
template<Expression E1, Expression E2, int Id, ConstantType Exponent, bool constants_compare_equal>
requires (is_constant_v<E1> && is_less_Sum_v<E2, Power<Symbol<Id>, Exponent>>)
struct is_less_Sum<Product<E1, E2>, Power<Symbol<Id>, Exponent>, constants_compare_equal> : /*std::true_type*/ ValueType<13> { };

// Compare Product with Product.
template<Expression E1, Expression E2, Expression E3, Expression E4, bool constants_compare_equal>
requires (
    (!is_constant_v<E1> && !is_constant_v<E3> &&
     (is_less_Sum_v<E1, E3> || (!is_less_Sum_v<E3, E1> && is_less_Sum_v<E2, E4>))) ||
    (!is_constant_v<E1> && is_constant_v<E3> && is_less_Sum_v<Product<E1, E2>, E4>) ||
    (is_constant_v<E1> && !is_constant_v<E3> && is_less_Sum_v<E2, Product<E3, E4>>) ||
    (is_constant_v<E1> && is_constant_v<E3> && is_less_Sum_v<E2, E4>))
struct is_less_Sum<Product<E1, E2>, Product<E3, E4>, constants_compare_equal> : /*std::true_type*/ ValueType<14> { };

template<Expression Base1, ConstantType Exponent1, Expression E2, bool constants_compare_equal>
requires (!is_constant_v<E2> && !is_symbol_v<E2> && !is_power_v<E2> && !is_product_v<E2> && !is_exponentiation_v<E2> && !is_multiplication_v<E2>)
struct is_less_Sum<Exponentiation<Base1, Exponent1>, E2, constants_compare_equal> : /*std::true_type*/ ValueType<15> { };

template<Expression Base1, ConstantType Exponent1, Expression E1, Expression E2, bool constants_compare_equal>
requires (is_constant_v<E1> && is_less_Sum_v<Exponentiation<Base1, Exponent1>, E2>)
struct is_less_Sum<Exponentiation<Base1, Exponent1>, Product<E1, E2>, constants_compare_equal> : /*std::true_type*/ ValueType<16> { };

template<Expression Base1, ConstantType Exponent1, Expression Base2, ConstantType Exponent2, bool constants_compare_equal>
requires (is_less_Sum_v<Base1, Base2> || (!is_less_Sum_v<Base2, Base1> && is_less_Sum<Exponent1, Exponent2, false>::value != 0))
struct is_less_Sum<Exponentiation<Base1, Exponent1>, Exponentiation<Base2, Exponent2>, constants_compare_equal> : /*std::true_type*/ ValueType<17> { };

template<Expression Base, ConstantType Exponent, Expression E1, Expression E2, bool constants_compare_equal>
requires (!is_constant_v<E1>)
struct is_less_Sum<Exponentiation<Base, Exponent>, Multiplication<E1, E2>, constants_compare_equal> : /*std::true_type*/ ValueType<18> { };

template<Expression E1, Expression E2, Expression E3, bool constants_compare_equal>
requires (!is_constant_v<E1> && !is_constant_v<E3> && !is_symbol_v<E3> && !is_power_v<E3> && !is_product_v<E3> && !is_exponentiation_v<E3> && !is_multiplication_v<E3>)
struct is_less_Sum<Multiplication<E1, E2>, E3, constants_compare_equal> : /*std::true_type*/ ValueType<19> { };

template<int Enumerator, int Denominator, Expression E2, Expression E3, bool constants_compare_equal>
requires (!is_multiplication_v<E3> && is_less_Sum_v<E2, E3>)
struct is_less_Sum<Multiplication<Constant<Enumerator, Denominator>, E2>, E3, constants_compare_equal> : /*std::true_type*/ ValueType<20> { };

template<Expression E1, int Enumerator, int Denominator, Expression E2, bool constants_compare_equal>
requires (!is_constant_v<E1> && !is_product_v<E1> && is_less_Sum_v<E1, E2>)
struct is_less_Sum<E1, Multiplication<Constant<Enumerator, Denominator>, E2>, constants_compare_equal> : /*std::true_type*/ ValueType<19> { };

template<Expression E1, Expression E2, Expression E3, Expression E4, bool constants_compare_equal>
requires(!is_constant_v<E3> &&
    ((!is_constant_v<E1> && (is_less_Sum_v<E1, E3> || (!is_less_Sum_v<E3, E1> && is_less_Sum_v<E2, E4>))) ||
     (is_constant_v<E1> && !is_constant_v<E3> && is_less_Sum_v<E2, Multiplication<E3, E4>>)))
struct is_less_Sum<Multiplication<E1, E2>, Multiplication<E3, E4>, constants_compare_equal> : /*std::true_type*/ ValueType<21> { };

template<Expression E1, Expression E2, bool constants_compare_equal>
requires (!is_constant_v<E2> && !is_symbol_v<E2> && !is_power_v<E2> && !is_product_v<E2> && !is_exponentiation_v<E2> && !is_multiplication_v<E2> && !is_sin_v<E2>)
struct is_less_Sum<Sin<E1>, E2, constants_compare_equal> : /*std::true_type*/ ValueType<22> { };

template<Expression E1, Expression E2, bool constants_compare_equal>
requires (is_less_Sum_v<E1, E2>)
struct is_less_Sum<Sin<E1>, Sin<E2>, constants_compare_equal> : /*std::true_type*/ ValueType<23> { };

template<Expression E1, Expression E2, bool constants_compare_equal>
requires (!is_constant_v<E2> && !is_symbol_v<E2> && !is_power_v<E2> && !is_product_v<E2> && !is_exponentiation_v<E2> && !is_multiplication_v<E2> && !is_sin_v<E2> && !is_cos_v<E2>)
struct is_less_Sum<Cos<E1>, E2, constants_compare_equal> : /*std::true_type*/ ValueType<22> { };

template<Expression E1, Expression E2, bool constants_compare_equal>
requires (is_less_Sum_v<E1, E2>)
struct is_less_Sum<Cos<E1>, Cos<E2>, constants_compare_equal> : /*std::true_type*/ ValueType<23> { };

// Compare two Sum's.
template<Expression E1, Expression E2, Expression E3, Expression E4, bool constants_compare_equal>
requires (is_less_Sum_v<E1, E3> || (!is_less_Sum_v<E3, E1> && is_less_Sum_v<E2, E4>))
struct is_less_Sum<Sum<E1, E2>, Sum<E3, E4>, constants_compare_equal> : /*std::true_type*/ ValueType<24> { };

} // namespace symbolic
