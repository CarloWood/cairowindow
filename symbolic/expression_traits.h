#pragma once

#include "Expression.h"
#include "is_product.h"         // constant, symbol, power and product.
#include "is_sum.h"
#include "is_exponentiation.h"
#include "is_multiplication.h"
#include "is_sin.h"
#include "is_cos.h"
#include "is_less_Sum.h"
#include "is_less_Product.h"
#include "is_less_Multiplication.h"

namespace symbolic {

template<Expression E>
struct get_exponent
{
  using type = Constant<1, 1>;
};

template<Expression E>
using get_exponent_t = typename get_exponent<E>::type;

template<Expression E>
struct get_base
{
  using type = E;
};

template<Expression E>
using get_base_t = typename get_base<E>::type;

template<Expression E>
struct get_constant_factor
{
  using type = Constant<1, 1>;
};

template<Expression E>
using get_constant_factor_t = typename get_constant_factor<E>::type;

template<int Enumerator1, int Denominator1>
struct get_constant_factor<Constant<Enumerator1, Denominator1>>
{
  using type = Constant<Enumerator1, Denominator1>;
};

template<int Enumerator1, int Denominator1, Expression E2>
struct get_constant_factor<Product<Constant<Enumerator1, Denominator1>, E2>>
{
  using type = Constant<Enumerator1, Denominator1>;
};

template<Expression E1, Expression E2>
struct get_constant_factor<Multiplication<E1, E2>>
{
  using type = get_constant_factor_t<E1>;
};

template<Expression E>
requires (!is_constant_v<E>)
struct get_nonconstant_factor
{
  using type = E;
};

template<Expression E>
using get_nonconstant_factor_t = typename get_nonconstant_factor<E>::type;

template<int Enumerator1, int Denominator1, Expression E2>
struct get_nonconstant_factor<Product<Constant<Enumerator1, Denominator1>, E2>>
{
  using type = E2;
};

} // namespace symbolic

#include "multiply.h"

namespace symbolic {

template<Expression E1, Expression E2>
struct get_nonconstant_factor<Multiplication<E1, E2>>
{
  static consteval auto eval()
  {
    if constexpr (is_constant_v<E1>)
      return get_nonconstant_factor_t<E2>::instance();
    else
      return multiply_t<get_nonconstant_factor_t<E1>, get_nonconstant_factor_t<E2>>::instance();
  }

 public:
  using type = decltype(eval());
};

} // namespace symbolic
