#pragma once

#include "IdRange.h"
#include "Constant.h"
#include "Symbol.h"
#include "Power.h"
#include "Expression.h"

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
requires (is_symbol_v<E1> || is_power_v<E1> || is_constant_v<E1>)
class Product : public ExpressionTag
{
 public:
  using arg1_type = E1;
  using arg2_type = E2;

  static constexpr precedence s_precedence = is_minus_one_v<E1> ? precedence::negation : precedence::product;
  static constexpr IdRange<E1::id_range.begin, E2::id_range.end> id_range{};

 private:
  E1 arg1_;
  E2 arg2_;

 public:
  inline constexpr Product(E1 const& arg1, E2 const& arg2);

  constexpr E1 const& arg1() const { return arg1_; }
  constexpr E2 const& arg2() const { return arg2_; }

#ifdef SYMBOLIC_PRINTING
  void print_on(std::ostream& os) const
  {
    if constexpr (s_precedence == precedence::negation)
    {
      bool need_parens = needs_parens(arg2_.s_precedence, s_precedence);
      os << "-";
      if (need_parens)
        os << '(';
      arg2_.print_on(os);
      if (need_parens)
        os << ')';
    }
    else
    {
      bool need_parens = needs_parens(arg1_.s_precedence, s_precedence, before);
      if (need_parens)
        os << '(';
      arg1_.print_on(os);
      if (need_parens)
        os << ')';
      os << " * ";
      need_parens = needs_parens(arg2_.s_precedence, s_precedence, after);
      if (need_parens)
        os << '(';
      arg2_.print_on(os);
      if (need_parens)
        os << ')';
    }
  }
#endif
};

template<typename T>
struct is_product : std::false_type { };

template<typename T>
constexpr bool is_product_v = is_product<T>::value;

template<Expression E1, Expression E2>
struct is_product<Product<E1, E2>> : std::true_type { };

template<typename E>
concept ProductType = is_product_v<E>;

template<Expression E1, Expression E2, Expression E3, Expression E4>
constexpr bool is_same_expression(Product<E1, E2> const& arg1, Product<E3, E4> const& arg2)
{
  return std::is_same_v<E1, E3> && std::is_same_v<E2, E4>;
}

template<Expression E1, Expression E2, Expression E3, Expression E4>
requires (is_less_v<E1, E3> || (!is_less_v<E3, E1> && is_less_v<E2, E4>))
struct is_less<Product<E1, E2>, Product<E3, E4>> : std::true_type { };

template<Expression E1, Expression E2, Expression E3>
requires (!is_constant_v<E3> && !is_symbol_v<E3> && !is_power_v<E3> && !is_product_v<E3>)
struct is_less<Product<E1, E2>, E3> : std::true_type { };

template<Expression E1, Expression E2>
struct is_product_less : std::false_type { };

template<Expression E1, Expression E2>
constexpr bool is_product_less_v = is_product_less<E1, E2>::value;

template<int Enumerator1, int Denominator1, Expression E2>
requires (!is_constant_v<E2>)
struct is_product_less<Constant<Enumerator1, Denominator1>, E2> : std::true_type { };

template<int Id1, SymbolPowerType E2>
requires (Id1 < E2::s_id)
struct is_product_less<Symbol<Id1>, E2> : std::true_type { };

template<SymbolType E1, int Enumerator1, int Denominator1, SymbolPowerType E2>
requires (E1::s_id < E2::s_id)
struct is_product_less<Power<E1, Enumerator1, Denominator1>, E2> : std::true_type { };

template<Expression E1, Expression E2>
requires (is_symbol_v<E1> || is_power_v<E1> || is_constant_v<E1>)
constexpr Product<E1, E2>::Product(E1 const& arg1, E2 const& arg2) : arg1_(arg1), arg2_(arg2)
{
  if constexpr (is_constant_v<E1>)
    static_assert(!E1::is_zero() && !E1::is_one(), "A Product should never contain zero or one as factor.");
  if constexpr (is_constant_v<E2>)
    static_assert(!E2::is_zero() && !E2::is_one(), "A Product should never contain zero or one as factor.");
  static_assert(is_symbol_v<E2> || is_power_v<E2> || is_product_v<E2>, "The second factor of a Product can only be a Symbol, Power or another Product.");
  if constexpr (!is_product_v<E2>)
    static_assert(is_product_less_v<E1, E2>, "The first argument of a Product must be less than the second argument.");
  else
    static_assert(is_product_less_v<E1, typename E2::arg1_type>,
        "The first argument of a Product must be less than the first argument of the second argument, if that is a Product.");
}

template<Expression E1, Expression E2>
constexpr auto multiply_unequals(E1 const& arg1, E2 const& arg2)
{
  if constexpr (is_product_v<E2>)
    static_assert(is_product_less_v<E1, typename E2::arg1_type>, "Only call multiply_unequals for E1 < E2.");
  else
    static_assert(is_product_less_v<E1, E2>, "Only call multiply_unequals for E1 < E2.");
  static_assert(!is_product_v<E1>, "Only call multiply_unequals with a non-product arg1.");

  if constexpr (is_constant_v<E1>)
  {
    if constexpr (E1::is_zero())
      return constant<0>();
    else if constexpr (E1::is_one())
      return arg2;
    else
      return Product{arg1, arg2};
  }
  else
    return Product{arg1, arg2};
}

// The design of the following multiplication operator is as follows:
//
// Prerequisite is that all Product have the form {non-Product, Expression}
// and that they are already ordered; in other words, each Product has
// a cannonical form: if they represent the same product, their types
// will be the same.
//
// This also means that every Product can be thought of as a chaining
// of non-Product's, because if in the above Expression is a Product,
// it itself begins again with a non-Product.
//
// Lets call the non-Products involved, a, b, c, d etc. Then a Product
// can be thought of as an ordered series: b * (c * (f * (h * x))).
//
// For Symbols the ordering is determined by their Id. There can only
// be one Constant because those contract eachother; a Constant is
// always "less" than a Symbol with respect to product ordering.
// Hence, in the above (only) 'b' can be a Constant, but that is not
// necessary.
//
// The ordering is determined by the smallest Symbol Id in the Expression
// (or -1 for a Constant). However, each Expression also references the
// largest Symbol/Constant Id (plus one) as `end`. This allows us to quickly
// determine if a product has an overlapping range or not.
//
// In the comments below I use '*' for an operator* operation and an 'x'
// for a Product operation. Thus (a x b x c) really is (a x (b x c))
// (or rather, Product{a, Product{b, c}}).
//
// Note that any of the symbols in the comments can be a Symbol, or a Power
// of that Symbol. Thus (a x b x c) can also be a^(3/2) x b^7 x c^(-1)).
//
// In order to be able to combine constants and powers of the same symbol,
// the ordering of factors must be as follows: constant < symbol/power,
// symbol/power(x) < symbol/power(y), etc.

template<Expression E1, Expression E2>
constexpr auto operator*(E1 const& arg1, E2 const& arg2)
{
  if constexpr (!is_product_v<E1> && !is_product_v<E2>)
  {
    if constexpr (is_product_less_v<E1, E2>)
      return multiply_unequals(arg1, arg2);
    else if constexpr (is_product_less_v<E2, E1>)
      return multiply_unequals(arg2, arg1);
    else
      return make_power(arg1, arg2);
  }
  else if constexpr (!is_product_v<E1>)
  {
    if constexpr (is_product_less_v<E1, typename E2::arg1_type>)
      return multiply_unequals(arg1, arg2);
    else if constexpr (is_product_less_v<typename E2::arg1_type, E1>)
      return multiply_unequals(arg2.arg1(), arg1 * arg2.arg2());
    else
      return make_power(arg1, arg2.arg1()) * arg2.arg2();
  }
  else if constexpr (!is_product_v<E2>)
  {
    if constexpr (is_product_less_v<typename E1::arg1_type, E2>)
      return multiply_unequals(arg1.arg1(), arg1.arg2() * arg2);
    else if constexpr (is_product_less_v<E2, typename E1::arg1_type>)
      return multiply_unequals(arg2, arg1);
    else
      return make_power(arg2, arg1.arg1()) * arg1.arg2();
  }
  // Both are a Product.
  else if constexpr (is_product_less_v<typename E1::arg1_type, typename E2::arg1_type>)
    return multiply_unequals(arg1.arg1(), arg1.arg2() * arg2);
  else if constexpr (is_product_less_v<typename E2::arg1_type, typename E1::arg1_type>)
    return multiply_unequals(arg2.arg1(), arg1 * arg2.arg2());
  else
    return make_power(arg1.arg1(), arg2.arg1()) * arg1.arg2() * arg2.arg2();
}

template<Expression E1, Expression E2>
constexpr auto inverse(Product<E1, E2> const& arg)
{
  auto arg1_inverse = inverse(arg.arg1());
  auto arg2_inverse = inverse(arg.arg2());
  if constexpr (is_less_v<decltype(arg1_inverse), decltype(arg2_inverse)>)
    return Product{arg1_inverse, arg2_inverse};
  else
    return Product{arg2_inverse, arg1_inverse};                                 // Note arg2_inverse can't be a product in this case,
                                                                                // because arg1_inverse isn't and !is_less_v<...>.
}

template<Expression E1, Expression E2>
requires (!is_constant_v<E1> || !is_constant_v<E2>)
constexpr auto operator/(E1 const& arg1, E2 const& arg2)
{
  return arg1 * inverse(arg2);
}

template<ProductType E1, ConstantType E2>
constexpr auto operator^(E1 const& arg1, E2 const& exponent)
{
  return Product{arg1.arg1()^exponent, arg1.arg2()^exponent};
}

// Negation.
template<Expression E>
constexpr auto operator-(E const& expression)
{
  return constant<-1>() * expression;
}

} // namespace symbolic
