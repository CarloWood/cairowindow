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
constexpr auto operator+(E1 const&, E2 const&)
{
  return add_t<E1, E2>::instance();
}

// Subtraction.
template<Expression E1, Expression E2>
constexpr auto operator-(E1 const&, E2 const&)
{
  return add_t<E1, negate_t<E2>>::instance();
}

// Exponentiation.
template<Expression Base, ConstantType Exponent>
constexpr auto operator^(Base const&, Exponent const&)
{
  return exponentiate_t<Base, Exponent>::instance();
}

} // namespace symbolic
