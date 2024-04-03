#pragma once

#include "Constant.h"
#include "Product.h"
#include "Multiplication.h"

namespace symbolic {

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

} // namespace symbolic
