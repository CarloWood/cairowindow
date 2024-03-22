#pragma once

#include "expression_traits.h"
#include "add.h"
#include "negate.h"
#include "precedence.h"
#include "IdRange.h"
#include <algorithm>

namespace symbolic {

template<Expression E1, Expression E2>
class Sum : public ExpressionTag
{
  static_assert(!is_sum_v<E1>, "The first term of a Sum must not be a Sum itself.");
  static_assert(!is_constant_v<E2>, "The second term of a Sum is not allowed to be a constant.");
  static_assert(!is_constant_zero_v<E1>, "Don't construct a Sum with a zero constant.");

 public:
  using arg1_type = E1;
  using arg2_type = E2;

  static constexpr precedence s_precedence = precedence::sum;
  static constexpr IdRange<std::min(E1::id_range.begin, E2::id_range.begin), std::max(E1::id_range.end, E2::id_range.end)> id_range{};

 public:
  static Sum instance() { return {}; }

#ifdef SYMBOLIC_PRINTING
  static void print_on(std::ostream& os);
#endif
};

template<Expression E1, Expression E2, Expression E3, Expression E4>
constexpr bool is_same_expression(Sum<E1, E2> const& arg1, Sum<E3, E4> const& arg2)
{
  return std::is_same_v<E1, E3> && std::is_same_v<E2, E4>;
}

#ifdef SYMBOLIC_PRINTING
//static
template<Expression E1, Expression E2>
void Sum<E1, E2>::print_on(std::ostream& os)
{
  bool need_parens = needs_parens(E1::s_precedence, s_precedence, before);
  if (need_parens)
    os << '(';
  E1::print_on(os);
  if (need_parens)
    os << ')';
  if constexpr (is_product_v<E2>)
  {
    if constexpr (is_constant_less_than_zero_v<typename E2::arg1_type>)
    {
      os << " - ";
      using Term = negate_t<E2>;
      need_parens = needs_parens(Term::s_precedence, precedence::difference, after);
      if (need_parens)
        os << '(';
      Term::print_on(os);
      if (need_parens)
        os << ')';
      return;
    }
  }
  else if constexpr (is_sum_v<E2>)
  {
    if constexpr (is_product_v<typename E2::arg1_type>)
    {
      if constexpr (is_constant_less_than_zero_v<typename E2::arg1_type::arg1_type>)
      {
        os << " - ";
        add_t<negate_t<typename E2::arg1_type>, typename E2::arg2_type>::print_on(os);
        return;
      }
    }
  }
  os << " + ";
  need_parens = needs_parens(E2::s_precedence, s_precedence, after);
  if (need_parens)
    os << '(';
  E2::print_on(os);
  if (need_parens)
    os << ')';
}
#endif

} // namespace symbolic
