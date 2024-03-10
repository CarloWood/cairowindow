#pragma once

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

 private:
  E const& expression_;

 public:
  consteval Negation(E const& expression) : expression_(expression) { }

  consteval E const& operator-() const { return expression_; }

#ifdef CWDEBUG
  bool needs_parens(before_or_after position, precedence prec) const { return prec < s_precedence; }
  bool needs_parens(precedence prec) const { return prec < s_precedence; }

  void print_on(std::ostream& os) const
  {
    bool need_parens = expression_.needs_parens(s_precedence);
    os << "-!";
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
consteval auto operator-(E const& expression)
{
  return Negation<E>{expression};
}

} // namespace symbolic
