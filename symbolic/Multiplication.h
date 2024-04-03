#pragma once

#include "IdRange.h"
#include "precedence.h"
#include "Expression.h"
#include "is_constant.h"

namespace symbolic {

template<Expression E1, Expression E2>
class Multiplication : public ExpressionTag
{
 public:
  using arg1_type = E1;
  using arg2_type = E2;

  static constexpr precedence s_precedence = is_constant_minus_one_v<E1> ? precedence::negation : precedence::product;
  static constexpr IdRange<E1::id_range.begin, E2::id_range.end> id_range{};

 public:
  static constexpr Multiplication instance() { return {}; }

#ifdef SYMBOLIC_PRINTING
  static void print_on(std::ostream& os)
  {
    if constexpr (s_precedence == precedence::negation)
    {
      bool need_parens = needs_parens(E2::s_precedence, s_precedence);
      os << "-";
      if (need_parens)
        os << '(';
      E2::print_on(os);
      if (need_parens)
        os << ')';
    }
    else
    {
      bool need_parens = needs_parens(E1::s_precedence, s_precedence, before);
      if (need_parens)
        os << '(';
      E1::print_on(os);
      if (need_parens)
        os << ')';
      os << " * ";
      need_parens = needs_parens(E2::s_precedence, s_precedence, after);
      if (need_parens)
        os << '(';
      E2::print_on(os);
      if (need_parens)
        os << ')';
    }
  }
#endif
};

} // namespace symbolic
