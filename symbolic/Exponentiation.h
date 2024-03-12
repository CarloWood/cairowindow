#pragma once

#include "Expression.h"
#include "Constant.h"

namespace symbolic {

template<Expression E1, int Enumerator, int Denominator>
class Exponentiation : public ExpressionTag
{
 public:
  static constexpr precedence s_precedence = precedence::exponentiation;
  static constexpr auto s_exponent = constant<Enumerator, Denominator>();

 private:
  E1 base_;

 public:
  constexpr Exponentiation(E1 const& base) : base_(base) { }

  consteval E1 const& base() const { return base_; }

#ifdef SYMBOLIC_PRINTING
  constexpr bool needs_parens(before_or_after position, precedence prec) const { return prec < s_precedence; }
  constexpr bool needs_parens(precedence prec) const { return prec < s_precedence; }

  void print_on(std::ostream& os) const
  {
    bool need_parens = base_.needs_parens(s_precedence);
    if (need_parens)
      os << '(';
    base_.print_on(os);
    if (need_parens)
      os << ')';
    os << '^';
    constexpr bool need_parens2 = s_exponent.needs_parens(s_precedence);
    if constexpr (need_parens2)
      os << '(';
    s_exponent.print_on(os);
    if constexpr (need_parens2)
      os << ')';
  }
#endif
};

} // namespace symbolic
