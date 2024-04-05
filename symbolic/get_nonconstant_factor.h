#pragma once

#include "multiply_fwd.h"
#include "Constant.h"
#include "Product.h"
#include "Multiplication.h"

namespace symbolic {

// As a special case - this will return the Constant when E is a Constant (currently used in is_less_Sum).
template<Expression E>
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

template<Expression E1, Expression E2>
struct get_nonconstant_factor<Multiplication<E1, E2>>
{
  static consteval auto eval()
  {
    static_assert(std::is_same_v<E2, get_nonconstant_factor_t<E2>>, "Huh - second argument of a Multiplication should never contain a constant factor!");
    if constexpr (is_constant_v<E1>)
      return get_nonconstant_factor_t<E2>::instance();
    else
      return multiply_t<get_nonconstant_factor_t<E1>, get_nonconstant_factor_t<E2>>::instance();
  }

 public:
  using type = decltype(eval());
};

} // namespace symbolic
