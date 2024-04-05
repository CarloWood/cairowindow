#pragma once

#include "Expression.h"
#include "precedence.h"
#include "is_constant.h"
#include "is_symbol.h"
#include "get_exponent.h"
#include "get_base.h"

namespace symbolic {

template<Expression Base, ConstantType Exponent>
requires (!is_symbol_v<Base> && !is_constant_zero_v<Exponent> && !is_constant_one_v<Exponent>)
class Exponentiation : public ExpressionTag
{
 public:
  using base_type = Base;
  using exponent_type = Exponent;

  // Used for is_less_Sum.
  using arg1_type = Base;
  using arg2_type = Exponent;

  static constexpr precedence s_precedence = precedence::power;
  static constexpr auto id_range = Base::id_range;

 public:
  static constexpr Exponentiation instance() { return {}; }

#ifdef SYMBOLIC_PRINTING
 public:
  static void print_on(std::ostream& os)
  {
    bool need_parens = needs_parens(Base::s_precedence, s_precedence, before);
    if (need_parens)
      os << '(';
    Base::print_on(os);
    if (need_parens)
      os << ')';
    os << '^';
    need_parens = needs_parens(Exponent::s_precedence, s_precedence, after);
    if (need_parens)
      os << '(';
    Exponent::print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<Expression Base, ConstantType Exponent>
struct get_exponent<Exponentiation<Base, Exponent>>
{
  using type = Exponent;
};

template<Expression Base, ConstantType Exponent>
struct get_base<Exponentiation<Base, Exponent>>
{
  using type = Base;
};

} // namespace symbolic
