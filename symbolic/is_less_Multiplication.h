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

template<typename T>
concept NonProductLevelType = !ProductLevelType<T>;

template<int Enumerator1, int Denominator1, NonProductLevelType E2>
struct is_less_Multiplication<Constant<Enumerator1, Denominator1>, E2> : std::true_type { };

template<int Id, NonProductLevelType E2>
struct is_less_Multiplication<Symbol<Id>, E2> : std::true_type { };

template<int Id, ConstantType Exponent, NonProductLevelType E2>
struct is_less_Multiplication<Power<Symbol<Id>, Exponent>, E2> : std::true_type { };

template<Expression E1, Expression E2, NonProductLevelType E3>
struct is_less_Multiplication<Product<E1, E2>, E3> : std::true_type { };

// Non-Product types are currently:
// - Sum                }
// - Multiplication     } -- SME
// - Exponentiation     }
// - Sin
// - Cos

template<Expression E1, Expression E2, NonProductLevelType E3>
requires (!is_sum_v<E3>)
struct is_less_Multiplication<Sum<E1, E2>, E3> : std::true_type { };

template<Expression E1, Expression E2, Expression E3, Expression E4>
requires (is_less_Multiplication_v<E1, E3> || (!is_less_Multiplication_v<E3, E1> && is_less_Multiplication_v<E2, E4>))
struct is_less_Multiplication<Sum<E1, E2>, Sum<E3, E4>> : std::true_type { };

template<Expression E1, Expression E2, Expression E3, Expression E4>
requires (is_less_Multiplication_v<E1, E3> || (!is_less_Multiplication_v<E3, E1> && is_less_Multiplication_v<E2, E4>))
struct is_less_Multiplication<Multiplication<E1, E2>, Multiplication<E3, E4>> : std::true_type { };

template<Expression E1, Expression E2, Expression Base, ConstantType Exponent>
requires (is_less_Multiplication_v<E1, Base>)
struct is_less_Multiplication<Multiplication<E1, E2>, Exponentiation<Base, Exponent>> : std::true_type { };

template<Expression Base, ConstantType Exponent, Expression E1, Expression E2>
requires (is_less_Multiplication_v<Base, E1>)
struct is_less_Multiplication<Exponentiation<Base, Exponent>, Multiplication<E1, E2>> : std::true_type { };

template<Expression Base1, ConstantType Exponent1, Expression Base2, ConstantType Exponent2>
requires (is_less_Multiplication_v<Base1, Base2>)
struct is_less_Multiplication<Exponentiation<Base1, Exponent1>, Exponentiation<Base2, Exponent2>> : std::true_type { };

template<typename T>
concept FunctionType = NonProductLevelType<T> && !is_sum_v<T> && !is_multiplication_v<T> && !is_exponentiation_v<T>;

template<Expression E1, Expression E2, FunctionType E3>
requires (is_less_Multiplication_v<E1, E3>)
struct is_less_Multiplication<Multiplication<E1, E2>, E3> : std::true_type { };

template<Expression Base, ConstantType Exponent, FunctionType E3>
requires (is_less_Multiplication_v<Base, E3>)
struct is_less_Multiplication<Exponentiation<Base, Exponent>, E3> : std::true_type { };

template<FunctionType E1, Expression E2, Expression E3>
requires (is_less_Multiplication_v<E1, E2>)
struct is_less_Multiplication<E1, Multiplication<E2, E3>> : std::true_type { };

template<FunctionType E1, Expression Base, ConstantType Exponent>
requires (is_less_Multiplication_v<E1, Base>)
struct is_less_Multiplication<E1, Exponentiation<Base, Exponent>> : std::true_type { };

template<Expression E1, Expression E2>
requires (is_less_Multiplication_v<E1, E2>)
struct is_less_Multiplication<Sin<E1>, Sin<E2>> : std::true_type { };

template<Expression E1, Expression E2>
requires (is_less_Multiplication_v<E1, E2>)
struct is_less_Multiplication<Cos<E1>, Cos<E2>> : std::true_type { };

template<Expression E1, Expression E2>
struct is_less_Multiplication<Sin<E1>, Cos<E2>> : std::true_type { };

} // namespace symbolic
