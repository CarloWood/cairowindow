#pragma once

#include "multiply.h"
#include "Constant.h"
#include <utility>

namespace symbolic {

template<Expression Base, ConstantType Exponent, IsProduct is_product =
  is_symbol_v<Base> || is_power_v<Base> || is_product_v<Base> ? isProduct : not_a_Product>
struct exponentiate;

template<Expression Base, ConstantType Exponent>
using exponentiate_t = typename exponentiate<Base, Exponent>::type;

template<Expression E, IsProduct is_product>
struct exponentiate<E, Constant<0, 1>, is_product>
{
  static_assert(!is_constant_zero_v<E>, "Zero to the power zero is undefined.");
  using type = Constant<1, 1>;
};

template<typename Exponent>
concept NonTrivialExponent = is_constant_v<Exponent> && !is_constant_zero_v<Exponent> && !is_constant_one_v<Exponent>;

template<ConstantType Base, NonTrivialExponent Exponent>
struct exponentiate<Base, Exponent, not_a_Product>
{
  static_assert(Exponent::s_denominator == 1, "Fractional exponents of constants are not (yet) supported.");

 private:
  static consteval auto eval()
  {
    if constexpr ((Exponent::s_enumerator & 1) == 1)
      return multiply_t<Base, exponentiate_t<Base, Constant<Exponent::s_enumerator - 1, 1>>>::instance();
    else
    {
      using Sqrt = exponentiate_t<Base, Constant<Exponent::s_enumerator / 2, 1>>;
      return multiply_t<Sqrt, Sqrt>::instance();
    }
  }

 public:
  using type = decltype(eval());
};

template<Expression E, IsProduct is_product>
struct exponentiate<E, Constant<1, 1>, is_product>
{
  using type = E;
};

template<int Id, NonTrivialExponent Exponent>
struct exponentiate<Symbol<Id>, Exponent, isProduct>
{
  using type = Power<Symbol<Id>, Exponent>;
};

template<SymbolType S, ConstantType Exponent1, ConstantType Exponent2>
struct exponentiate<Power<S, Exponent1>, Exponent2, isProduct>
{
  using type = exponentiate_t<S, multiply_t<Exponent1, Exponent2>>;
};

template<Expression E1, Expression E2, NonTrivialExponent Exponent>
struct exponentiate<Product<E1, E2>, Exponent, isProduct>
{
  using type = Product<exponentiate_t<E1, Exponent>, exponentiate_t<E2, Exponent>>;
};

template<Expression E, ConstantType Exponent1, NonTrivialExponent Exponent2>
struct exponentiate<Exponentiation<E, Exponent1>, Exponent2, not_a_Product>
{
  using type = exponentiate_t<E, multiply_t<Exponent1, Exponent2>>;
};

template<Expression E, NonTrivialExponent Exponent>
requires (!is_constant_v<E> && !is_exponentiation_v<E>)
struct exponentiate<E, Exponent, not_a_Product>
{
  using type = Exponentiation<E, Exponent>;
};

} // namespace symbolic
