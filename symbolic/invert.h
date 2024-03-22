#pragma once

#include "negate.h"

namespace symbolic {

// Forward declarations.

template<Expression Base, ConstantType Exponent>
class Exponentiation;

// Define invert for all types.

template<Expression E>
struct invert;

template<Expression E>
using invert_t = typename invert<E>::type;

template<int Enumerator, int Denominator>
requires (Enumerator != 0)
struct invert<Constant<Enumerator, Denominator>>
{
  static constexpr int sign = Enumerator > 0 ? 1 : -1;
  using type = Constant<sign * Denominator, sign * Enumerator>;
};

template<int Id>
struct invert<Symbol<Id>>
{
  using type = Power<Symbol<Id>, Constant<-1, 1>>;
};

template<Expression Base, ConstantType Exponent>
struct invert<Power<Base, Exponent>>
{
  using type =
    /*if*/ std::conditional_t<is_constant_minus_one_v<Exponent>,
      Base,
    /*else*/
      Power<Base, negate_t<Exponent>> >;
};

template<Expression E1, Expression E2>
struct invert<Product<E1, E2>>
{
  using arg1_inverse_t = invert_t<E1>;
  using arg2_inverse_t = invert_t<E2>;

  using type =
    /*if*/ std::conditional_t<is_less_Product<arg1_inverse_t, arg2_inverse_t>::value,
      Product<arg1_inverse_t, arg2_inverse_t>,
    /*else*/
      Product<arg2_inverse_t, arg1_inverse_t> >;        // Note arg2_inverse_t can't be a Product in this case,
                                                        // because arg1_inverse_t isn't and !is_less_Product_v<...>.
};

template<Expression E1, Expression E2>
struct invert<Sum<E1, E2>>
{
  using type = Exponentiation<Sum<E1, E2>, Constant<-1, 1>>;
};

template<Expression Base, ConstantType Exponent>
struct invert<Exponentiation<Base, Exponent>>
{
  using type =
  /*if*/ std::conditional_t<is_constant_minus_one_v<Exponent>,
    Base,
  /*else*/
    Exponentiation<Base, negate_t<Exponent>> >;
};

template<Expression E>
auto old_invert(E const&)
{
  return invert_t<E>::instance();
}

} // namespace symbolic
