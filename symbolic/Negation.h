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
  using base_type = E;

  static constexpr precedence s_precedence = precedence::negation;
  static constexpr auto id_range = E::id_range;

 private:
  E expression_;

 public:
  inline constexpr Negation(E const& expression);

  constexpr E const& operator-() const { return expression_; }

#ifdef SYMBOLIC_PRINTING
  void print_on(std::ostream& os) const
  {
    bool need_parens = needs_parens(expression_.s_precedence, s_precedence);
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
requires (!ConstantType<E>)
constexpr Negation<E>::Negation(E const& expression) : expression_(expression)
{
  if constexpr (is_product_v<E>)
    static_assert(!is_constant_v<typename E::arg1_type>, "Should never get a Negation of a Product of a Constant times something.");
  static_assert(!is_negation_v<E>, "It should never happen that a Negation of a Negation is constructed!");
}

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
