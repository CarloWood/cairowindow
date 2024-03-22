#pragma once

#include "Symbol.h"
#include "Power.h"
#include "Product.h"
#include "Sum.h"
#include "Exponentiation.h"
#include "negate.h"
#include "multiply.h"
#include "invert.h"
#include "add.h"
#include "exponentiate.h"

namespace symbolic {

// Multiplication.
template<Expression E1, Expression E2>
constexpr auto operator*(E1 const& arg1, E2 const& arg2)
{
  return multiply_t<E1, E2>::instance();
}

// Division.
template<Expression E1, Expression E2>
requires (!is_constant_v<E1> || !is_constant_v<E2>)
constexpr auto operator/(E1 const& arg1, E2 const& arg2)
{
  return multiply_t<E1, invert_t<E2>>::instance();
}

// Negation.
template<Expression E>
constexpr auto operator-(E const&)
{
  return negate_t<E>::instance();
}

// Addition.
template<Expression E1, Expression E2>
auto operator+(E1 const&, E2 const&)
{
  return add_t<E1, E2>::instance();
}

// Subtraction.
template<Expression E1, Expression E2>
auto operator-(E1 const&, E2 const&)
{
  return add_t<E1, negate_t<E2>>::instance();
}

// Exponentiation.
template<SymbolType E1, ConstantType E2>
constexpr auto operator^(E1 const& base, E2 const& exponent)
{
  return Power<E1, E2>::instance();
}

template<PowerType Base, ConstantType Exponent>
constexpr auto operator^(Base const& power, Exponent const& exponent)
{
  return Power<typename Base::base_type, multiply_t<typename Base::exponent_type, Exponent>>::instance();
}

template<Expression E1, Expression E2, ConstantType Exponent>
constexpr auto operator^(Product<E1, E2> const&, Exponent const&)
{
  return Product<exponentiate_t<E1, Exponent>, exponentiate_t<E2, Exponent>>::instance();
}
template<SumType Base, ConstantType Exponent>

auto operator^(Base const& sum, Exponent const& exponent)
{
  return Exponentiation<Base, Exponent>{sum};
}

template<ExponentiationType Base, ConstantType Exponent>
auto operator^(Base const& exponentiation, Exponent const& exponent)
{
  return Exponentiation<typename Base::base_type, multiply_t<typename Base::exponent_type, Exponent>>::instance();
}

} // namespace symbolic
