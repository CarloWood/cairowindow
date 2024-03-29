#pragma once

#include "multiply_fwd.h"
#include "precedence.h"
#include "IdRange.h"
#include <numeric>              // std::gcd
#ifdef SYMBOLIC_PRINTING
#include "utils/has_print_on.h"
#include <iostream>
#endif
#include "debug.h"

namespace symbolic {

#ifdef SYMBOLIC_PRINTING
using utils::has_print_on::operator<<;
#endif

// Forward declarations (needed for the friend declarations).

template<int E, int D>
consteval auto constant();

template<int E, int D, int exponent>
consteval auto operator^(Constant<E, D>, Constant<exponent, 1>);

template<int E, int D>
consteval auto operator-(Constant<E, D>);

// A rational constant.
//
template<int Enumerator, int Denominator>
class Constant : public ExpressionTag
{
  static_assert(Denominator > 0, "The denominator should always be positive.");
  static_assert(std::gcd(Enumerator, Denominator) == 1, "A Constant must be constructed in canonical form.");

 public:
  static constexpr precedence s_precedence =
    Denominator != 1 ? precedence::product : Enumerator < 0 ? precedence::negation : precedence::constant;
  static constexpr int s_enumerator = Enumerator;
  static constexpr int s_denominator = Denominator;
  static constexpr IdRange<-1, 0> id_range{};

 private:
  // Use this free function to create Constant objects.
  template<int E, int D>
  friend consteval auto constant();

  // Negation.
  template<int E, int D>
  friend consteval auto operator-(Constant<E, D>);

  template<int E, int D, int exponent>
  friend consteval auto operator^(Constant<E, D>, Constant<exponent, 1>);

 public:
  static constexpr Constant instance() { return {}; }

#ifdef SYMBOLIC_PRINTING
 public:
  static void print_on(std::ostream& os)
  {
    os << Enumerator;
    if constexpr (Denominator > 1)
      os << "/" << Denominator;
  }
#endif
};

} // namespace symbolic

#include "is_constant.h"
#include "canonical_constant.h"

namespace symbolic {

template<Expression E1, Expression E2>
class add;

template<int Enumerator1, int Denominator1, int Enumerator2, int Denominator2>
struct add<Constant<Enumerator1, Denominator1>, Constant<Enumerator2, Denominator2>>
{
  using type = canonical_constant_t<Enumerator1 * Denominator2 + Enumerator2 * Denominator1, Denominator1 * Denominator2>;
};

template<int Enumerator1, int Denominator1, int Enumerator2, int Denominator2>
struct multiply<Constant<Enumerator1, Denominator1>, Constant<Enumerator2, Denominator2>, not_a_Product>
{
  using type = canonical_constant_t<Enumerator1 * Enumerator2, Denominator1 * Denominator2>;
};

// Create a rational constant.
//
// Usage:
//
// constexpr auto my_constant = symbolic::constant<13, 37>();
//
// The returned type is canonicalized: the denominator is always > 0
// and the GCD of enumerator and denominator will be 1.
//
template<int E, int D = 1>
consteval auto constant()
{
  return canonical_constant_t<E, D>::instance();
}

// Negation.
template<int E, int D>
consteval auto operator-(Constant<E, D>)
{
  return Constant<-E, D>{};
}

// Addition.
template<int E1, int D1, int E2, int D2>
consteval auto operator+(Constant<E1, D1>, Constant<E2, D2>)
{
  return constant<E1 * D2 + E2 * D1, D1 * D2>();
}

// Subtraction.
template<int E1, int D1, int E2, int D2>
consteval auto operator-(Constant<E1, D1>, Constant<E2, D2>)
{
  return constant<E1 * D2 - E2 * D1, D1 * D2>();
}

// Multiplication.
template<int E1, int D1, int E2, int D2>
consteval auto operator*(Constant<E1, D1>, Constant<E2, D2>)
{
  return constant<E1 * E2, D1 * D2>();
}

// Division.
template<int E1, int D1, int E2, int D2>
requires (E2 != 0)
consteval auto operator/(Constant<E1, D1>, Constant<E2, D2>)
{
  return constant<E1 * D2, D1 * E2>();
}

// Exponentiation with an integer.
template<int E, int D, int exponent>
consteval auto operator^(Constant<E, D>, Constant<exponent, 1>)
{
  if constexpr (exponent == 0)
  {
    static_assert(E != 0, "0^0 is undefined");
    return Constant<1, 1>{};
  }
  else if constexpr (exponent == 1)
    return Constant<E, D>{};
  else if constexpr (exponent < 0)
  {
    static_assert(E != 0, "division by zero is undefined");
    return Constant<D, E>{}^Constant<-exponent, 1>{};
  }
  else if constexpr ((exponent & 1) == 1)
    return Constant<E, D>{} * (Constant<E, D>{}^Constant<exponent - 1, 1>{});
  else
    return constant<E * E, D * D>() ^ Constant<exponent / 2, 1>{};;
}

} // namespace symbolic
