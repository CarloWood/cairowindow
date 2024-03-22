#pragma once

//#include "Constant.h"
#include "expression_traits.h"
#include "precedence.h"

namespace symbolic {

template<Expression Base, ConstantType Exponent>
class Power : public ExpressionTag
{
  static_assert(is_symbol_v<Base>, "Base must be a Symbol.");

 public:
  using base_type = Base;
  using exponent_type = Exponent;
  static constexpr int s_id = Base::s_id;

  static constexpr precedence s_precedence = precedence::power;
  static constexpr auto id_range = Base::id_range;

  static Power instance() { return {}; }

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

template<SymbolType Base, ConstantType Exponent>
struct get_base<Power<Base, Exponent>>
{
  using type = Base;
};

} // namespace symbolic
