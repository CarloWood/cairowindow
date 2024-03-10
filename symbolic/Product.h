#pragma once

#include "Expression.h"

namespace symbolic {

template<Expression E1, Expression E2>
class Product : public ExpressionTag
{
 public:
  static constexpr precedence s_precedence = precedence::product;

 private:
  E1 const& arg1_;
  E2 const& arg2_;

 public:
  consteval Product(E1 const& arg1, E2 const& arg2) : arg1_(arg1), arg2_(arg2) { }

  consteval E1 const& arg1() const { return arg1_; }
  consteval E2 const& arg2() const { return arg2_; }

#ifdef CWDEBUG
  bool needs_parens(before_or_after position, precedence prec) const { return prec < s_precedence; }
  bool needs_parens(precedence prec) const { return prec < s_precedence; }

  void print_on(std::ostream& os) const
  {
    bool need_parens = arg1_.needs_parens(before, s_precedence);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " * ";
    need_parens = arg2_.needs_parens(after, s_precedence);
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<Expression E1, Expression E2>
consteval auto operator*(E1 const& arg1, E2 const& arg2)
{
  static constexpr expression_order_less<E2, E1> eol(arg2, arg1);
  static constexpr bool foo = eol();
  if constexpr (foo)
    return Product{arg2, arg1};
  else
    return Product{arg1, arg2};
}

} // namespace symbolic
