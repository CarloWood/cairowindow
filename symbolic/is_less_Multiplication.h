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

template<Expression E1, Expression E2, bool is_product = ProductLevelType<E1> && ProductLevelType<E2>>
struct is_less_Multiplication : std::false_type { };

template<Expression E1, Expression E2>
constexpr bool is_less_Multiplication_v = is_less_Multiplication<E1, E2>::value;

// If the arguments to is_less_Multiplication are both Product-level types, then use is_less_Product.
template<ProductLevelType E1, ProductLevelType E2>
struct is_less_Multiplication<E1, E2, true> : is_less_Product<E1, E2> { };

template<Expression E1, Expression E2>
requires ((!is_exponentiation_v<E1> && !is_exponentiation_v<E2>) &&
          ((is_multiplication_v<E2> && !is_multiplication_v<E1>) || is_less_exact_v<E1, E2>))
struct is_less_Multiplication<E1, E2, false> : TrueType<10> { };

template<Expression Base, ConstantType Exponent, Expression E2>
requires (!is_exponentiation_v<E2> && is_less_exact_v<Base, get_base_t<E2>>)
struct is_less_Multiplication<Exponentiation<Base, Exponent>, E2, false> : TrueType<11> { };

template<Expression E1, Expression Base, ConstantType Exponent>
requires (is_less_exact_v<get_base_t<E1>, Base>)
struct is_less_Multiplication<E1, Exponentiation<Base, Exponent>, false> : TrueType<12> { };

} // namespace symbolic
