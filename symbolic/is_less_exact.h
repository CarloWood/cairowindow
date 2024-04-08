#pragma once

#include "is_product.h"
#include "is_exponentiation.h"
#include "is_multiplication.h"
#include "is_sum.h"
#include "is_sin.h"
#include "is_cos.h"
#include "is_log.h"
#include "is_atan.h"

namespace symbolic {

// Two expressions are only "equal" when they differ at most a constant factor.
// Any other ordering is therefore arbitrary - as long as terms that only
// differ a constant factor are brought together.
//
// Neverless, we use the following ordering:
//
//   Constant < Symbol < Power < Product < Exponentiation < Multiplication < Sin < Cos < Log < Sum.
//
// where we can assume that a Power and Exponentiation never have an exponent of 0 or 1
// and the Product nor Multiplication begin with a constant (if they do, then the non-constant
// factor has to be compared instead).
//
// Note that Product < Multiplication in all cases, because it is impossible
// for the Multiplication to exist of a Constant times something that is less
// than a Product.
//
// In case two arguments are the same, the following holds:
//
//   1) two Constant's : compared by value (for 'exact' comparisons).
//   2) two Symbol's : compared by Id.
//   3) two Power's : first compared by Base then by exponent.
//   4) two Product's : first compared by first argument, then by second argument.
//   5) two Exponentiation's : first compared by Base then by exponent.
//   6) two Multiplication's : first compared by first argument, then by second argument.
//   7) two Sin : compared by argument.
//   8) two Cos : compared by argument.
//   9) two Log : compared by argument.
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
struct expression_order
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
    else if constexpr (is_log_v<E>)
      return unary_op | 4;                      // 01100
    else if constexpr (is_atan_v<E>)
      return unary_op | 5;                      // 01101
    // Binary operators:
    else if constexpr (is_sum_v<E>)
      return binary_op | 16;                    // 10010
    else
      static_assert(DependentFalse<E>::value, "Not implemented");
  }

  static constexpr int value = eval();
};

template<Expression E>
constexpr int expression_order_v = expression_order<E>::value;

enum class OpType
{
  special, binary, unary
};

template<Expression E>
constexpr bool is_binary_op = (expression_order_v<E> & binary_op);

template<Expression E>
constexpr bool is_unary_op = (expression_order_v<E> & unary_op);

template<Expression E1, Expression E2, OpType = (is_binary_op<E1> ? OpType::binary : is_unary_op<E1> ? OpType::unary : OpType::special)>
struct is_less_same_kind_exact : TrueType<0> { };

template<Expression E1, Expression E2>
constexpr bool is_less_same_kind_exact_v = is_less_same_kind_exact<E1, E2>::value;

template<Expression E1, Expression E2>
constexpr bool is_less_exact_v =
    (expression_order_v<E1> < expression_order_v<E2> || (expression_order_v<E1> == expression_order_v<E2> && is_less_same_kind_exact_v<E1, E2>));

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
requires (is_less_exact_v<typename E1::arg1_type, typename E2::arg1_type> ||
          (!is_less_exact_v<typename E2::arg1_type, typename E1::arg1_type> &&
           is_less_exact_v<typename E1::arg2_type, typename E2::arg2_type>))
struct is_less_same_kind_exact<E1, E2, OpType::binary> : TrueType<4> { };

// All unary operators are compared by their argument.
template<Expression E1, Expression E2>
requires (is_less_exact_v<typename E1::arg_type, typename E2::arg_type>)
struct is_less_same_kind_exact<E1, E2, OpType::unary> : TrueType<5> { };

} // namespace symbolic
