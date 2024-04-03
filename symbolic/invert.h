#pragma once

#include "negate.h"
#include "is_less_Product.h"

namespace symbolic {

// Forward declarations.

template<Expression Base, ConstantType Exponent>
requires (!is_constant_zero_v<Exponent> && !is_constant_one_v<Exponent>)
class Exponentiation;

// Define invert for all types.

template<Expression E>
struct invert
{
  using type = Exponentiation<E, Constant<-1, 1>>;
};

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
 private:
  static consteval auto eval()
  {
    if constexpr (is_constant_minus_one_v<Exponent>)
      return Base::instance();
    else
      return Power<Base, negate_t<Exponent>>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2>
struct invert<Product<E1, E2>>
{
  using arg1_inverse_t = invert_t<E1>;
  using arg2_inverse_t = invert_t<E2>;

 private:
  static consteval auto eval()
  {
    if constexpr (is_less_Product_v<arg1_inverse_t, arg2_inverse_t>)
      return Product<arg1_inverse_t, arg2_inverse_t>::instance();
    else
      return Product<arg2_inverse_t, arg1_inverse_t>::instance();       // Note arg2_inverse_t can't be a Product in this case,
                                                                        // because arg1_inverse_t isn't and !is_less_Product_v<...>.
  }

 public:
  using type = decltype(eval());
};

template<Expression Base, ConstantType Exponent>
struct invert<Exponentiation<Base, Exponent>>
{
 private:
  static consteval auto eval()
  {
    if constexpr (is_constant_minus_one_v<Exponent>)
     return Base::instance();
    else
     return Exponentiation<Base, negate_t<Exponent>>::instance();
  }

 public:
  using type = decltype(eval());
};

} // namespace symbolic
