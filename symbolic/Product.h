#pragma once

#include "IdRange.h"
#include "precedence.h"
#include "is_product.h"

namespace symbolic {

// A Product is a product of- optionally exponentiated- symbols and an optional constant;
// it has a canonical form which can be defined as:
//
//   ❮product❯                    ::= ❮product{-1,l}❯ | ❮product{f,l}❯ where f and l are integer identifiers and 0 <= f <= l.
//   ❮product{-1,l}❯              ::= Product<❮constant❯, ❮power{l}❯> | Product<❮constant❯, ❮constant-free-product{f,l}❯>
//   ❮product{f,l}❯               ::= ❮constant-free-product{f,l}❯
//   ❮power{s}❯                   ::= ❮symbol{s}❯ | Power<❮symbol{s}❯, ❮constant-not-zero-or-one❯> where 0 <= s.
//   ❮constant❯                   ::= Constant<E, D> (in its canonical form: gcd(E, D) = 1 and D > 0).
//   ❮constant-not-zero-or-one❯   ::= ❮constant❯ - Constant<0, 1> - Constant<1, 1>
//   ❮constant-free-product{f,l}❯ ::= Product<❮power{f}❯, ❮power{l}❯> | Product<❮power{f}❯, ❮constant-free-product{g,l}❯> where 0 <= f < g < l.
//
// Examples:
//
//   -x                                 Product<Constant<-1, 1>, Symbol<0>>
//   2 * x^2                            Product<Constant<2, 1>, Power<Symbol<0>, 2, 1>>
//   x * y                              Product<Symbol<0>, Symbol<1>>
//   -3/2 * x^(-5/4) * y^(7/3)          Product<Constant<-3, 2>, Product<Power<Symbol<0>, -5, 4>, Power<Symbol<1>, 7, 3>>>
//
// Note that each symbol only occurs once and the product is sorted by Symbol Id.
//
template<Expression E1, Expression E2>
requires (((is_constant_v<E1> && !is_constant_zero_v<E1> && !is_constant_one_v<E1>) ||
           is_symbol_v<E1> || is_power_v<E1>) &&
          (is_symbol_v<E2> || is_power_v<E2> || is_product_v<E2>))
class Product : public ExpressionTag
{
 public:
  using arg1_type = E1;
  using arg2_type = E2;

  static constexpr precedence s_precedence = is_constant_minus_one_v<E1> ? precedence::negation : precedence::product;
  static constexpr IdRange<E1::id_range.begin, E2::id_range.end> id_range{};

 public:
  inline constexpr Product();

  static constexpr Product instance() { return {}; }

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

template<Expression E1, Expression E2, Expression E3, Expression E4>
constexpr bool is_same_expression(Product<E1, E2> const& arg1, Product<E3, E4> const& arg2)
{
  return std::is_same_v<E1, E3> && std::is_same_v<E2, E4>;
}

} // namespace symbolic

#include "is_less_Product.h"

namespace symbolic {

template<Expression E1, Expression E2>
requires (((is_constant_v<E1> && !is_constant_zero_v<E1> && !is_constant_one_v<E1>) ||
           is_symbol_v<E1> || is_power_v<E1>) &&
          (is_symbol_v<E2> || is_power_v<E2> || is_product_v<E2>))
constexpr Product<E1, E2>::Product()
{
  if constexpr (!is_product_v<E2>)
    static_assert(is_less_Product_v<E1, E2>, "The first argument of a Product must be less than the second argument.");
  else
    static_assert(is_less_Product_v<E1, typename E2::arg1_type>,
        "The first argument of a Product must be less than the first argument of the second argument, if that is a Product.");
}

} // namespace symbolic
