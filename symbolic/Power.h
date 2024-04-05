#pragma once

#include "precedence.h"
#include "is_symbol.h"
#include "is_constant.h"
#include "get_exponent.h"

namespace symbolic {

template<SymbolType Base, ConstantType Exponent>
requires (!is_constant_zero_v<Exponent> && !is_constant_one_v<Exponent>)
class Power : public ExpressionTag
{
 public:
  using base_type = Base;
  using exponent_type = Exponent;

  // Used for is_less_Sum.
  using arg1_type = Base;
  using arg2_type = Exponent;

  static constexpr int s_id = Base::s_id;

  static constexpr precedence s_precedence = precedence::power;
  static constexpr auto id_range = Base::id_range;

  static constexpr Power instance() { return {}; }

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

template<SymbolType Base, ConstantType Exponent>
struct get_exponent<Power<Base, Exponent>>
{
  using type = Exponent;
};

} // namespace symbolic
