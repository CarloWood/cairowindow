#pragma once

#include "expression_traits.h"

#ifndef MULTIPLY_H
#error "multiply<> should be fully defined before using it."
#endif

namespace symbolic {

template<Expression E, SymbolType S>
struct derivative
{
  using type = Constant<0, 1>;
};

template<Expression E, SymbolType S>
using derivative_t = typename derivative<E, S>::type;

template<int Id>
struct derivative<Symbol<Id>, Symbol<Id>>
{
  using type = Constant<1, 1>;
};

template<int Id, ConstantType Exponent>
struct derivative<Power<Symbol<Id>, Exponent>, Symbol<Id>>
{
  using type = multiply_t<Exponent, exponentiate_t<Symbol<Id>, add_t<Exponent, Constant<-1, 1>>>>;
};

template<Expression E1, Expression E2, int Id>
requires (Product<E1, E2>::id_range.begin <= Id && Id < Product<E1, E2>::id_range.end)
struct derivative<Product<E1, E2>, Symbol<Id>>
{
  using type = add_t<multiply_t<E1, derivative_t<E2, Symbol<Id>>>, multiply_t<derivative_t<E1, Symbol<Id>>, E2>>;
};

template<Expression E1, Expression E2, int Id>
struct derivative<Sum<E1, E2>, Symbol<Id>>
{
  using type = add_t<derivative_t<E1, Symbol<Id>>, derivative_t<E2, Symbol<Id>>>;
};

template<Expression E, ConstantType Exponent, int Id>
struct derivative<Exponentiation<E, Exponent>, Symbol<Id>>
{
  using type = multiply_t<derivative_t<E, Symbol<Id>>, multiply_t<Exponent, exponentiate_t<E, add_t<Exponent, Constant<-1, 1>>>>>;
};

template<Expression T, SymbolType S>
auto differentiate(T const&, S const&)
{
  return derivative_t<T, S>::instance();
}


} // namespace symbolic
