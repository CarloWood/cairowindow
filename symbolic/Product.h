#pragma once

#include "Expression.h"
#include "Exponentiation.h"

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
  constexpr bool needs_parens(before_or_after position, precedence prec) const { return prec < s_precedence; }
  constexpr bool needs_parens(precedence prec) const { return prec < s_precedence; }

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
  if constexpr (expression_order_less(std::type_identity<E1>{}, std::type_identity<E2>{}))
    return Product{arg1, arg2};
  else if constexpr (expression_order_less(std::type_identity<E2>{}, std::type_identity<E1>{}))
    return Product{arg2, arg1};
  else
  {
    ASSERT(is_same_expression(arg1, arg2));
    return Exponentiation<E1, 2, 1>{arg1};
  }
}

} // namespace symbolic
