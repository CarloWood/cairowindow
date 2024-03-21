#pragma once

#include "Expression.h"

namespace symbolic {

template<Expression E, int Enumerator, int Denominator>
class Exponentiation : public ExpressionTag
{
 public:
  using base_type = E;
  static constexpr auto s_exponent = constant<Enumerator, Denominator>();

  static constexpr precedence s_precedence = precedence::power;
  static constexpr auto id_range = E::id_range;

 private:
  E base_;

 public:
  Exponentiation(E const& base) : base_(base)
  {
    static_assert(!is_symbol_v<E>, "Use Power to exponentiate symbols.");
  }

  E const& base() const { return base_; }

#ifdef SYMBOLIC_PRINTING
  void print_on(std::ostream& os) const
  {
    bool need_parens = needs_parens(base_.s_precedence, s_precedence, before);
    if (need_parens)
      os << '(';
    base_.print_on(os);
    if (need_parens)
      os << ')';
    os << '^';
    need_parens = needs_parens(s_exponent.s_precedence, s_precedence, after);
    if (need_parens)
      os << '(';
    s_exponent.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<typename T>
struct is_exponentiation : std::false_type { };

template<typename T>
constexpr bool is_exponentiation_v = is_exponentiation<T>::value;

template<Expression E, int Enumerator, int Denominator>
struct is_exponentiation<Exponentiation<E, Enumerator, Denominator>> : std::true_type { };

template<typename E>
concept ExponentiationType = is_exponentiation_v<E>;

template<ExponentiationType E>
consteval auto get_exponent()
{
  return E::s_exponent;
}

template<SumType E1, ConstantType E2>
auto operator^(E1 const& sum, E2 exponent)
{
  return Exponentiation<E1, E2::s_enumerator, E2::s_denominator>{sum};
}

template<ExponentiationType E1, ConstantType E2>
auto operator^(E1 const& exponentiation, E2 exponent)
{
  constexpr auto new_exponent = constant<E1::s_exponent.s_enumerator * E2::s_enumerator, E1::s_exponent.s_denominator * E2::s_denominator>();
  return Exponentiation<typename E1::base_type, new_exponent.s_enumerator, new_exponent.s_denominator>{exponentiation.base()};
}

template<Expression E1, int Enumerator, int Denominator>
constexpr auto inverse(Exponentiation<E1, Enumerator, Denominator> const& exponentiation)
{
  if constexpr (Enumerator == -1 && Denominator == 1)
    return exponentiation.base();
  else
    return Exponentiation<E1, -Enumerator, Denominator>{exponentiation.base()};
}

} // namespace symbolic
