#pragma once

#include "Expression.h"
#include "precedence.h"
#include "IdRange.h"
#include <numeric>
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
template<int Enumerator, int Denominator>
class Constant;

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

  // Constructor.
  consteval Constant() = default;

 public:
  static consteval bool is_less_than_zero() { return Enumerator < 0; }
  static consteval bool is_minus_one() { return Enumerator == -1 && Denominator == 1; }
  static consteval bool is_zero() { return Enumerator == 0; }
  static consteval bool is_one() { return Enumerator == 1 && Denominator == 1; }

#ifdef SYMBOLIC_PRINTING
 public:
  void print_on(std::ostream& os) const
  {
    os << Enumerator;
    if constexpr (Denominator > 1)
      os << "/" << Denominator;
  }
#endif
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
  constexpr int abs_enumerator = E < 0 ? -E : E;
  constexpr int abs_denominator = D < 0 ? -D : D;
  constexpr int gcd = std::gcd(abs_enumerator, abs_denominator);
  if constexpr (D < 0)
    return Constant<-E / gcd, -D / gcd>{};
  else
    return Constant<E / gcd, D / gcd>{};
}

static consteval auto make_one() { return constant<1, 1>(); }
static consteval auto make_minus_one() { return constant<-1, 1>(); }

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

template<typename T>
struct is_constant : std::false_type {};

template<typename T>
constexpr bool is_constant_v = is_constant<T>::value;

template<int Enumerator, int Denominator>
struct is_constant<Constant<Enumerator, Denominator>> : std::true_type {};

template<typename E>
concept ConstantType = is_constant_v<E>;

template<int Enumerator, int Denominator>
constexpr auto inverse(Constant<Enumerator, Denominator> const&)
{
  return constant<Denominator, Enumerator>();
}

} // namespace symbolic
