#pragma once

#include "Product.h"
#include "Expression.h"
#include "Constant.h"
#include "precedence.h"
#include "utils/derived_from_template.h"

namespace symbolic {

template<Expression E>
requires (!ConstantType<E>)
class Negation : public ExpressionTag
{
 public:
  static constexpr precedence s_precedence = precedence::negation;
  static constexpr auto id_range = E::id_range;

 private:
  E expression_;

 public:
  constexpr Negation(E const& expression) : expression_(expression)
  {
    if constexpr (is_product_v<E>)
    {
      static_assert(!is_constant_v<std::decay_t<decltype(expression.arg1())>>, "This should never happen!");
    }
  }

  constexpr E const& operator-() const { return expression_; }

#ifdef SYMBOLIC_PRINTING
  constexpr bool needs_parens(before_or_after position, precedence prec) const { return prec < s_precedence; }
  constexpr bool needs_parens(precedence prec) const { return prec < s_precedence; }

  void print_on(std::ostream& os) const
  {
    bool need_parens = expression_.needs_parens(s_precedence);
    os << "-";
    if (need_parens)
      os << '(';
    expression_.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<typename T>
struct is_negation : std::false_type {};

template<typename T>
constexpr bool is_negation_v = is_negation<T>::value;

template<Expression E>
struct is_negation<Negation<E>> : std::true_type {};

template<typename E>
concept NegationType = is_negation_v<E>;

template<Expression E>
requires (!NegationType<E>)
constexpr auto operator-(E const& expression)
{
  if constexpr (is_product_v<E>)
  {
    if constexpr (is_constant_v<typename E::arg1_type>)
      return Product{constant<-E::arg1_type::s_enumerator, E::arg1_type::s_denominator>(), expression.arg2()};
    else
      return Negation<E>{expression};
  }
  else
    return Negation<E>{expression};
}

template<Expression E1, Expression E2>
requires (!is_negation_v<E2>)
constexpr auto operator*(Negation<E1> const& arg1, E2 const& arg2)
{
  return -(-arg1 * arg2);
}

// Push negations out of products to the front, or combine them with a Constant.
template<Expression E1, Expression E2>
constexpr auto operator*(E1 const& arg1, Negation<E2> const& arg2)
{
  return -arg1 * -arg2;
}

template<Expression E>
constexpr auto inverse(Negation<E> const& arg)
{
  return Negation{inverse(-arg)};
}

} // namespace symbolic
