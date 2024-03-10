#pragma once

#include "Expression.h"
#include "precedence.h"
#include <numeric>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include <iostream>
#endif
#include "debug.h"

namespace symbolic {
using utils::has_print_on::operator<<;

// A rational constant.
//
template<int Enumerator, int Denominator>
class Constant : public ExpressionTag
{
 public:
  static constexpr precedence s_precedence = Denominator != 1 ? precedence::division : Enumerator < 0 ? precedence::negation : precedence::constant;

 private:
  // Use this free function to create Constant objects.
  template<int enumerator, int denominator>
  friend constexpr auto constant();

  // Negation.
  template<int E, int D>
  friend auto operator-(Constant<E, D>);

  template<int E, int D, int exponent>
  friend auto operator^(Constant<E, D>, Constant<exponent, 1>);

  // Constructor.
  constexpr Constant() = default;

#ifdef CWDEBUG
 public:
  bool needs_parens(before_or_after position, precedence prec) const { return prec < s_precedence; }
  bool needs_parens(precedence prec) const { return prec < s_precedence; }

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
template<int enumerator, int denominator>
constexpr auto constant()
{
  constexpr int abs_enumerator = enumerator < 0 ? -enumerator : enumerator;
  constexpr int abs_denominator = denominator < 0 ? -denominator : denominator;
  constexpr int gcd = std::gcd(abs_enumerator, abs_denominator);
  if constexpr (denominator < 0)
    return Constant<-enumerator / gcd, -denominator / gcd>{};
  else
    return Constant<enumerator / gcd, denominator / gcd>{};
}

// Negation.
template<int E, int D>
auto operator-(Constant<E, D>)
{
  return Constant<-E, D>{};
}

// Addition.
template<int E1, int D1, int E2, int D2>
auto operator+(Constant<E1, D1>, Constant<E2, D2>)
{
  DoutEntering(dc::notice, "operator+(Constant<" << E1 << ", " << D1 << ">, Constant<" << E2 << ", " << D2 << ">)");
  return constant<E1 * D2 + E2 * D1, D1 * D2>();
}

// Subtraction.
template<int E1, int D1, int E2, int D2>
auto operator-(Constant<E1, D1>, Constant<E2, D2>)
{
  DoutEntering(dc::notice, "operator-(Constant<" << E1 << ", " << D1 << ">, Constant<" << E2 << ", " << D2 << ">)");
  return constant<E1 * D2 - E2 * D1, D1 * D2>();
}

// Multiplication.
template<int E1, int D1, int E2, int D2>
auto operator*(Constant<E1, D1>, Constant<E2, D2>)
{
  DoutEntering(dc::notice, "operator*(Constant<" << E1 << ", " << D1 << ">, Constant<" << E2 << ", " << D2 << ">)");
  return constant<E1 * E2, D1 * D2>();
}

// Division.
template<int E1, int D1, int E2, int D2>
requires (E2 != 0)
auto operator/(Constant<E1, D1>, Constant<E2, D2>)
{
  DoutEntering(dc::notice, "operator/(Constant<" << E1 << ", " << D1 << ">, Constant<" << E2 << ", " << D2 << ">)");
  return constant<E1 * D2, D1 * E2>();
}

// Exponentiation with an integer.
template<int E, int D, int exponent>
auto operator^(Constant<E, D>, Constant<exponent, 1>)
{
  DoutEntering(dc::notice, "operator^(Constant<" << E << ", " << D << ">, Constant<" << exponent << ", 1>)");
  if constexpr (exponent == 0)
  {
    static_assert(E != 0, "0^0 is undefined");
    return Constant<1, 1>{};
  }
  if constexpr (exponent == 1)
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

} // namespace symbolic
