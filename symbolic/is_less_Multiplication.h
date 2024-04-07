#pragma once

#include "is_less_Product.h"
#include "is_multiplication.h"
#include "is_product.h"
#include "is_exponentiation.h"
#include "is_sum.h"
#include "Constant.h"
#include "Symbol.h"
#include "Power.h"
#include "Product.h"
#include "Sum.h"
#include "Multiplication.h"
#include "Exponentiation.h"
#include "functions.h"

namespace symbolic {

// Used to compare factors without exponentiation shortcircuiting.
// This is the same as is_less_exact_v, but anything not a multiplication compares less than a multiplication.
template<Expression E1, Expression E2>
constexpr bool is_less_Multiplication_exact_v = (is_multiplication_v<E2> && !is_multiplication_v<E1>) || is_less_exact_v<E1, E2>;

template<Expression E1, Expression E2, bool is_product = ProductLevelType<E1> && ProductLevelType<E2>>
struct is_less_Multiplication : std::false_type { };

template<Expression E1, Expression E2>
constexpr bool is_less_Multiplication_v = is_less_Multiplication<E1, E2>::value;

// If the arguments to is_less_Multiplication are both Product-level types, then use is_less_Product.
template<ProductLevelType E1, ProductLevelType E2>
struct is_less_Multiplication<E1, E2, true> : is_less_Product<E1, E2> { };

// Otherwise (if at least one argument is not a Product-level type) - if neither are exponentiations
// (which are handled below) - then if the second argument is a multiplication and the first one isn't,
// the first one is considered "less" than the second one (this is required for the ordering of the
// arguments of a Multiplication, where the first argument is not allowed to be a Multiplication), or
// when the second argument is not a Multiplication, or the first one is, we use is_less_exact_v.
//
// This gives the following comparison table:
//
//   E2:             E1: Product       Multiplication
//   Product              -               false
//   Multiplication      true           is_less_exact_v<E1, E2>
//
// E1 * E2 * E3 * E4 * ...
//
// Multiplication<E1, Multiplication<E2, Multiplication<E3, Multiplication<E4, ...>>>>
// where E1, E2, E3, E4, ... are not Multiplication's.
//
// Furthermore, get_base_t<E1> < get_base_t<E2> < get_base_t<E3> < get_base_t<E4>, which
// is required in order to be able to combine factors with the same "base".
// For example, say get_base_t<E2> == get_base_t<E3>, then the type would have been converted into
//
// Multiplication<E1, Multiplication<Exponentiation<get_base_t<E2>, get_exponent_v<E2> + get_exponent_v<E3>>, Multiplication<E4, ...>>>
//
template<Expression E1, Expression E2>
requires ((!is_exponentiation_v<E1> && !is_exponentiation_v<E2>) && is_less_Multiplication_exact_v<E1, E2>)
struct is_less_Multiplication<E1, E2, false> : TrueType<10> { };

template<Expression Base, ConstantType Exponent, Expression E2>
requires (!is_exponentiation_v<E2> && is_less_Multiplication_exact_v<Base, get_base_t<E2>>)
struct is_less_Multiplication<Exponentiation<Base, Exponent>, E2, false> : TrueType<11> { };

template<Expression E1, Expression Base, ConstantType Exponent>
requires (is_less_Multiplication_exact_v<get_base_t<E1>, Base>)
struct is_less_Multiplication<E1, Exponentiation<Base, Exponent>, false> : TrueType<12> { };

} // namespace symbolic
