#pragma once

#include "is_less_exact.h"
#include "get_nonconstant_factor.h"

namespace symbolic {

// struct is_less_Sum<E1, E2>
//
// This struct determines if expression E1 should be sorted left of expression E2
// when adding them - where possible constant factors of E1 and E2 are ignored.
//
// This sorting is required to get terms that are the SAME next to eachother;
// for example, let E1 = 3 * X and E2 = 5 * X, then both will be sorted next
// to eachother such that they can be combined as in: 3 * X + 5 * X = (3 + 5) * X.
template<Expression E1, Expression E2, bool = is_constant_factor_free_v<E1> && is_constant_factor_free_v<E2> && is_less_exact_v<E1, E2>>
struct is_less_Sum : TrueType<0> { };

// Ignore leading constant factors.
template<Expression E1, Expression E2>
constexpr bool is_less_Sum_v = is_less_Sum<E1, E2>::value;

// If both compared expressions are constant-free (not a Product or Multiplication that contains a constant factor),
// then we might as well compare them using is_less_exact_v<> (the value of which we have as third parameter).
// The exception is when both expressions are constants; we must never return true in that case.
template<Expression E1, Expression E2>
requires (!(is_constant_v<E1> && is_constant_v<E2>))
struct is_less_Sum<E1, E2, true> : TrueType<6> { };

// If one of the expression is a Product that contains a constant factor, then
// compare the non-constant factor (ignoring the constant factor, or assuming it is 1).
template<int e, int d, Expression E1, Expression E2>
requires (is_constant_factor_free_v<E2> && is_less_exact_v<E1, E2>)
struct is_less_Sum<Product<Constant<e, d>, E1>, E2, false> : TrueType<6> { };

template<Expression E1, int e, int d, Expression E2>
requires (is_less_exact_v<get_nonconstant_factor_t<E1>, E2>)
struct is_less_Sum<E1, Product<Constant<e, d>, E2>, false> : TrueType<7> { };

// Same for a Multiplication that contains a constant factor.
//
// Note that in the case of a Multiplication, we can have a Constant factor
// that is *part* of the first argument, as part of a Product: Multiplication<Product<Constant<...>, ...>, ...>.

template<int e, int d, Expression E1, Expression E2>
requires (is_constant_factor_free_v<E2> && is_less_exact_v<E1, E2>)
struct is_less_Sum<Multiplication<Constant<e, d>, E1>, E2, false> : TrueType<8> { };

template<Expression E1, int e, int d, Expression E2>
requires (is_less_exact_v<get_nonconstant_factor_t<E1>, E2>)
struct is_less_Sum<E1, Multiplication<Constant<e, d>, E2>, false> : TrueType<9> { };

template<int e, int d, Expression E1, Expression E2, Expression E3>
requires (is_constant_factor_free_v<E2> && is_less_exact_v<get_nonconstant_factor_t<Multiplication<Product<Constant<e, d>, E3>, E1>>, E2>)
struct is_less_Sum<Multiplication<Product<Constant<e, d>, E3>, E1>, E2, false> : TrueType<10> { };

template<Expression E1, int e, int d, Expression E2, Expression E3>
requires (is_less_exact_v<get_nonconstant_factor_t<E1>, get_nonconstant_factor_t<Multiplication<Product<Constant<e, d>, E3>, E2>>>)
struct is_less_Sum<E1, Multiplication<Product<Constant<e, d>, E3>, E2>, false> : TrueType<11> { };

} // namespace symbolic
