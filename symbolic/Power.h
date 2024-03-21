#pragma once

#include "Expression.h"
#include "Constant.h"
#include <type_traits>

namespace symbolic {

template<SymbolType S, int Enumerator, int Denominator>
class Power : public ExpressionTag
{
 public:
  using base_type = S;
  static constexpr auto s_exponent = constant<Enumerator, Denominator>();
  static constexpr int s_id = S::s_id;

  static constexpr precedence s_precedence = precedence::power;
  static constexpr auto id_range = S::id_range;

 private:
  S base_;

 public:
  constexpr Power(S const& base) : base_(base) { }

  S const& base() const { return base_; }

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
struct is_power : std::false_type { };

template<typename T>
constexpr bool is_power_v = is_power<T>::value;

template<SymbolType S, int Enumerator, int Denominator>
struct is_power<Power<S, Enumerator, Denominator>> : std::true_type { };

template<typename S>
concept PowerType = is_power_v<S>;

template<typename E>
concept SymbolPowerType = is_symbol_v<E> || is_power_v<E>;

template<SymbolType S>
consteval Constant<1, 1> get_exponent()
{
  return constant<1, 1>();
}

template<PowerType P>
consteval auto get_exponent()
{
  return P::s_exponent;
}

template<SymbolPowerType E1, SymbolPowerType E2>
constexpr auto make_power(E1 const& arg1, E2 const& arg2)
{
  static_assert(E1::id_range.begin == E2::id_range.begin, "Can only combine powers of the same symbol!");
  constexpr auto total_exponent = get_exponent<E1>() + get_exponent<E2>();

  if constexpr (total_exponent.s_enumerator == 0)
    return constant<1, 1>();
  else if constexpr (total_exponent.s_enumerator == 1 && total_exponent.s_denominator == 1)
  {
    if constexpr (is_symbol_v<E1>)
      return arg1;
    else
      return arg1.base();
  }
  else if constexpr (is_symbol_v<E1>)
    return Power<E1, total_exponent.s_enumerator, total_exponent.s_denominator>{arg1};
  else
    return Power<typename E1::base_type, total_exponent.s_enumerator, total_exponent.s_denominator>{arg1.base()};
}

template<int Enumerator1, int Denominator1, int Enumerator2, int Denominator2>
constexpr auto make_power(Constant<Enumerator1, Denominator1> const&, Constant<Enumerator2, Denominator2> const&)
{
  return constant<Enumerator1 * Enumerator2, Denominator1 * Denominator2>();
}

template<SymbolType E1, ConstantType E2>
constexpr auto operator^(E1 const& base, E2 const& exponent)
{
  return Power<E1, E2::s_enumerator, E2::s_denominator>{base};
}

template<PowerType E1, ConstantType E2>
constexpr auto operator^(E1 const& power, E2 const& exponent)
{
  constexpr auto new_exponent = constant<E1::s_exponent.s_enumerator * E2::s_enumerator, E1::s_exponent.s_denominator * E2::s_denominator>();
  return Power<typename E1::base_type, new_exponent.s_enumerator, new_exponent.s_denominator>{power.base()};
}

template<int Id>
constexpr auto inverse(Symbol<Id> symbol)
{
  return Power<Symbol<Id>, -1, 1>{symbol};
}

template<SymbolType E1, int Enumerator, int Denominator>
constexpr auto inverse(Power<E1, Enumerator, Denominator> const& power)
{
  if constexpr (Enumerator == -1 && Denominator == 1)
    return power.base();
  else
    return Power<E1, -Enumerator, Denominator>{power.base()};
}

} // namespace symbolic
