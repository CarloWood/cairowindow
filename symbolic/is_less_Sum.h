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
#include "Exponentiation.h"
#include "Multiplication.h"
#include "Sum.h"
#include "functions.h"
#include "get_nonconstant_factor.h"

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
// and the Product nor Multiplication begin with a constant (if they do, then the non-constant
// factor has to be compared instead).
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
// Note that in the case of a Multiplication, we can have a Constant factor
// that is *part* of the first argument, as part of a Product: Multiplication<Product<Constant<...>, ...>, ...>.
//
#if 0
// This is used to detect which specialization was used.
template<int v>
struct TrueType
{
  static constexpr int value = v;
  consteval operator bool() { return v != 0; }
};
#else
template<int v>
struct TrueType : std::true_type { };

template<>
struct TrueType<0> : std::false_type { };
#endif

static constexpr int binary_op = 2;             //    v
static constexpr int unary_op = 8;              //  v

// A helper struct to determine the global order of available expression classes.
template<Expression E>
struct sum_order
{
  static consteval int eval()
  {
    // Special cases:
    if constexpr (is_constant_v<E>)
      return 0;                                 // 00000
    else if constexpr (is_symbol_v<E>)
      return 1;                                 // 00001
    // Binary operators:
    else if constexpr (is_power_v<E>)
      return binary_op;                         // 00010
    else if constexpr (is_product_v<E>)
      return binary_op | 1;                     // 00011
    else if constexpr (is_exponentiation_v<E>)
      return binary_op | 4;                     // 00110
    else if constexpr (is_multiplication_v<E>)
      return binary_op | 5;                     // 00111
    // Unary operators:
    else if constexpr (is_sin_v<E>)
      return unary_op;                          // 01000
    else if constexpr (is_cos_v<E>)
      return unary_op | 1;                      // 01001
    // Binary operators:
    else if constexpr (is_sum_v<E>)
      return binary_op | 16;                    // 10010
    else
      static_assert(DependentFalse<E>::value, "Not implemented");
  }

  static constexpr int value = eval();
};

template<Expression E>
constexpr int sum_order_v = sum_order<E>::value;

enum class OpType
{
  special, binary, unary
};

template<Expression E>
constexpr bool is_binary_op = (sum_order_v<E> & binary_op);

template<Expression E>
constexpr bool is_unary_op = (sum_order_v<E> & unary_op);

template<Expression E1, Expression E2, OpType = (is_binary_op<E1> ? OpType::binary : is_unary_op<E1> ? OpType::unary : OpType::special)>
struct is_less_same_kind_exact : TrueType<0> { };

template<Expression E1, Expression E2>
constexpr bool is_less_same_kind_exact_v = is_less_same_kind_exact<E1, E2>::value;

// bool is_less_Sum_exact_v<E1, E2>
//
// This variable determines if Foo<E1> should be sorted left of Foo<E2>,
// where possible constant factors of E1 and E2 are NOT ignored.
// Foo here is either a unary operator, or a binary operator were
// E1 is the first argument, or the first argument is equal and E1
// is the second argument.
//
// This sorting is required so that Foo<E1> and Foo<E2> are NOT
// considered equal when they can not be combined as above; for example
// 3 * Foo<E1> + 5 * Foo<E2> = (3 + 5) * Foo<E1> only when E1 == E2
// exactly (do not differ a constant factor), which is the case when
// Foo<C * E> != C * Foo<E>, where C is a constant not equal to 0 or 1.
//
// The remaining part of this comment shows why this hold for every
// expression class that we use:
//
// Lets consider all unary operators first:
//
// Symbol<C * Id> != C * Symbol<Id> (those would be different symbols).
// Sin<C * E> != C * Sin<E>
// Cos<C * E> != C * Cos<E>
//
// Binary operators:
//
// Power<C * B, Ex> != C * Power<B, Ex>
// Power<B, C * Ex> != C * Power<B, Ex>
// Exponentiation<C * B, Ex> != C * Exponentiation<B, Ex>
// Exponentiation<B, C * Ex> != C * Exponentiation<B, Ex>
// Sum<C * X, Y> != C * Sum<X, Y>       - unless Y is zero, which is not allowed for a Sum.
// Sum<X, C * Y> != C * Sum<X, Y>       - unless X is zero, which is not allowed for a Sum.
// Note that Sum<C * X, C * Y> == C * Sum<X, Y>, but a Sum existing of two expressions both
// with the same constant factor must still result in Sum<C * X, C * Y> (multiplications
// are distributed into Sums).
//
// However (let Q be a Symbol or Power, and E a Symbol, Power or Product),
//
// Product<C * Q, E> == C * Product<Q, E>       - but Product<C * Q, E> is not allowed (the first argument can not be a Product itself).
// Product<Q, C * E> == C * Product<Q, E>       - but Product<Q, C * E> is not allowed (that would be Product<C, Product<Q, E>>).
//
// Let M be a non-product type (Exponentiation or higher in order)
//
// Multiplication<C * Q, M> = C * Multiplication<Q, M>
// Multiplication<Q, C * M> = C * Multiplication<Q, M>  - but Multiplication<Q, C * M> is not allowed
//                                                        (that would be Multiplication<C, Multiplication<Q, E>>).
// In other words, while obviously
//
// 3 * Multiplication<E1, M> + 5 * Multiplication<E1, M> = (3 + 5) * Multiplication<E1, M>
//
// we can also simplify when E1 and E2 differ a constant factor (here E is the non-constant factor of E1 and E2 which is consider equal):
//
// 3 * Multiplication<C * E, M> + 5 * Multiplication<D * E, M> = (3 * C + 5 * D) * Multiplication<E, M>.
//
// However, this can't happen because 3 * Multiplication<C * E, M> = Multiplication<3, Multiplication<C * E, M>>
// and the second argument of a Multiplication must be constant factor free (aka, C and D will always be 1).
//

// Do not ignore leading constant factors.
//
// For an exact compare, we sort expressions using their `sum_order` (as determined by struct sum_order above)
// and use is_less_same_kind_exact_v if the two expression have the same sum_order value.
template<Expression E1, Expression E2>
constexpr bool is_less_Sum_exact_v =
    (sum_order_v<E1> < sum_order_v<E2> || (sum_order_v<E1> == sum_order_v<E2> && is_less_same_kind_exact_v<E1, E2>));

// Two constants are compared by their value.
template<int e1, int d1, int e2, int d2>
requires (e1 * d2 < e2 * d1)
struct is_less_same_kind_exact<Constant<e1, d1>, Constant<e2, d2>, OpType::special> : TrueType<2> { };

// Two symbols are compared by their id.
template<int id1, int id2>
requires (id1 < id2)
struct is_less_same_kind_exact<Symbol<id1>, Symbol<id2>, OpType::special> : TrueType<3> { };

// All binary operators are first compared by their first argument, and if that is equal by their second argument.
template<Expression E1, Expression E2>
requires (is_less_Sum_exact_v<typename E1::arg1_type, typename E2::arg1_type> ||
          (!is_less_Sum_exact_v<typename E2::arg1_type, typename E1::arg1_type> &&
           is_less_Sum_exact_v<typename E1::arg2_type, typename E2::arg2_type>))
struct is_less_same_kind_exact<E1, E2, OpType::binary> : TrueType<4> { };

// All unary operators are compared by their argument.
template<Expression E1, Expression E2>
requires (is_less_Sum_exact_v<typename E1::arg_type, typename E2::arg_type>)
struct is_less_same_kind_exact<E1, E2, OpType::unary> : TrueType<5> { };

// struct is_less_Sum<E1, E2>
//
// This struct determines if expression E1 should be sorted left of expression E2
// when adding them - where possible constant factors of E1 and E2 are ignored.
//
// This sorting is required to get terms that are the SAME next to eachother;
// for example, let E1 = 3 * X and E2 = 5 * X, then both will be sorted next
// to eachother such that they can be combined as in: 3 * X + 5 * X = (3 + 5) * X.
template<Expression E1, Expression E2, bool = is_constant_factor_free_v<E1> && is_constant_factor_free_v<E2> && is_less_Sum_exact_v<E1, E2>>
struct is_less_Sum : TrueType<0> { };

// Ignore leading constant factors.
template<Expression E1, Expression E2>
constexpr bool is_less_Sum_v = is_less_Sum<E1, E2>::value;

// If both compared expressions are constant-free (not a Product or Multiplication that contains a constant factor),
// then we might as well compare them using is_less_Sum_exact_v<> (the value of which we have as third parameter).
// The exception is when both expressions are constants; we must never return true in that case.
template<Expression E1, Expression E2>
requires (!(is_constant_v<E1> && is_constant_v<E2>))
struct is_less_Sum<E1, E2, true> : TrueType<6> { };

// If one of the expression is a Product that contains a constant factor, then
// compare the non-constant factor (ignoring the constant factor, or assuming it is 1).
template<int e, int d, Expression E1, Expression E2>
requires (is_constant_factor_free_v<E2> && is_less_Sum_exact_v<E1, E2>)
struct is_less_Sum<Product<Constant<e, d>, E1>, E2, false> : TrueType<6> { };

template<Expression E1, int e, int d, Expression E2>
requires (is_less_Sum_exact_v<get_nonconstant_factor_t<E1>, E2>)
struct is_less_Sum<E1, Product<Constant<e, d>, E2>, false> : TrueType<7> { };

// Same for a Multiplication that contains a constant factor.
template<int e, int d, Expression E1, Expression E2>
requires (is_constant_factor_free_v<E2> && is_less_Sum_exact_v<E1, E2>)
struct is_less_Sum<Multiplication<Constant<e, d>, E1>, E2, false> : TrueType<8> { };

template<Expression E1, int e, int d, Expression E2>
requires (is_less_Sum_exact_v<get_nonconstant_factor_t<E1>, E2>)
struct is_less_Sum<E1, Multiplication<Constant<e, d>, E2>, false> : TrueType<9> { };

template<int e, int d, Expression E1, Expression E2, Expression E3>
requires (is_constant_factor_free_v<E2> && is_less_Sum_exact_v<get_nonconstant_factor_t<Multiplication<Product<Constant<e, d>, E3>, E1>>, E2>)
struct is_less_Sum<Multiplication<Product<Constant<e, d>, E3>, E1>, E2, false> : TrueType<10> { };

template<Expression E1, int e, int d, Expression E2, Expression E3>
requires (is_less_Sum_exact_v<get_nonconstant_factor_t<E1>, get_nonconstant_factor_t<Multiplication<Product<Constant<e, d>, E3>, E2>>>)
struct is_less_Sum<E1, Multiplication<Product<Constant<e, d>, E3>, E2>, false> : TrueType<11> { };

} // namespace symbolic
