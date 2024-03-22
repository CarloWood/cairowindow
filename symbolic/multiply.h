#pragma once

#include "expression_traits.h"
#include "add.h"
#include "Product.h"

namespace symbolic {

template<Expression E1, Expression E2>
struct make_power;

template<Expression E1, Expression E2>
using make_power_t = typename make_power<E1, E2>::type;

template<Expression E1, Expression E2>
struct multiply_unequals
{
  static consteval bool sanity_check()
  {
    if constexpr (is_product_v<E2>)
      return is_less_Product_v<E1, typename E2::arg1_type>;
    else
      return is_less_Product_v<E1, E2>;
  }

  static_assert(!is_product_v<E1>, "Only call multiply_unequals with a non-product arg1.");
  static_assert(sanity_check(), "Only call multiply_unequals for E1 < E2.");

  using type =
    /*if*/ std::conditional_t<is_constant_zero_v<E1>,
      Constant<0, 1>,
    /*else if*/ std::conditional_t<is_constant_one_v<E1>,
      E2,
    /*else*/
      Product<E1, E2> >>;
};

template<Expression E1, Expression E2>
using multiply_unequals_t = typename multiply_unequals<E1, E2>::type;

// The design of the multiply specialization with isProduct set is as follows:
//
// Prerequisite is that all Product types have the form: Product<non-Product, Expression>,
// where non-Product can be a Constant, a Symbol or a Power (of a Symbol),
// and Expression can be that including another Product - but only one
// that doesn't begin with a Constant.
// Also the Products are already ordered; in other words, each Product has
// a cannonical form: if they represent the same product, their types
// will be the same.
//
// This also means that every Product can be thought of as a chaining
// of Symbols and/or Powers, prepended with an optional Constant, because if
// in the above Expression is a Product, it itself begins with a Symbol
// or Power. Each such factor is or a different Symbol or a Power of a symbol
// that are ordered by Id.
//
// For example, let the non-Products involved, a, b, c, d etc. Then a
// Product can be thought of as the ordered product: b * (c * (f * (h * x))).
//
// For Symbols the ordering is determined by their Id. There can only
// be one Constant because those contract each other; a Constant is
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
// for a Product operation. Thus (a x b x c) really is Product{a, Product{b, c}}.
//
// Note that any of the symbols in the comments can be a Symbol, or a Power
// of that Symbol. Thus (a x b x c) can also be a^(3/2) x b^7 x c^(-1)).
//
// In order to be able to combine constants and powers of the same symbol,
// the ordering of factors must be as follows: constant < symbol/power,
// symbol/power(x) < symbol/power(y), etc.

// This default is for when E1 nor E2 is a Product (and therefore a Symbol,
// Power or (one of them) a Constant). If one or both are a Product then
// one of the specializations, with isProduct, below will be used.
//
template<Expression E1, Expression E2>
struct multiply<E1, E2, isProduct>
{
 private:
  static consteval auto eval()
  {
    if constexpr (is_less_Product_v<E1, E2>)
      return multiply_unequals_t<E1, E2>::instance();
    else if constexpr (is_less_Product_v<E2, E1>)
      return multiply_unequals_t<E2, E1>::instance();
    else
      return make_power_t<E1, E2>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2, Expression E3>
struct multiply<E1, Product<E2, E3>, isProduct>
{
 private:
  static consteval auto eval()
  {
    if constexpr (is_less_Product_v<E1, E2>)
      return multiply_unequals_t<E1, Product<E2, E3>>::instance();
    else if constexpr (is_less_Product_v<E2, E1>)
      return multiply_unequals_t<E2, multiply_t<E1, E3>>::instance();
    else
      return multiply_t<make_power_t<E1, E2>, E3>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2, Expression E3>
struct multiply<Product<E1, E2>, E3, isProduct>
{
 private:
  static consteval auto eval()
  {
    if constexpr (is_less_Product_v<E1, E3>)
      return multiply_unequals_t<E1, multiply_t<E2, E3>>::instance();
    else if constexpr (is_less_Product_v<E3, E1>)
      return multiply_unequals_t<E3, Product<E1, E2>>::instance();
    else
      return multiply_t<make_power_t<E3, E1>, E2>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2, Expression E3, Expression E4>
struct multiply<Product<E1, E2>, Product<E3, E4>, isProduct>
{
  static constexpr bool E1_less_than_E3 = is_less_Product_v<E1, E3>;
  static constexpr bool E3_less_than_E1 = is_less_Product_v<E3, E1>;
  static constexpr bool E1_equal_E3 = !E1_less_than_E3 && !E3_less_than_E1;

 private:
  static consteval auto eval()
  {
    if constexpr (is_less_Product_v<E1, E3>)
      return multiply_unequals_t<E1, multiply_t<E2, Product<E3, E4>>>::instance();
    else if constexpr (is_less_Product_v<E3, E1>)
      return multiply_unequals_t<E3, multiply_t<Product<E1, E2>, E4>>::instance();
    else
      return multiply_t<multiply_t<make_power_t<E1, E3>, E2>, E4>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2, Expression E3>
struct multiply<E1, Sum<E2, E3>, not_a_Product>
{
  using type =
    add_t<multiply_t<E1, E2>, multiply_t<E1, E3>>;
};

template<Expression E1, Expression E2, Expression E3>
struct multiply<Sum<E1, E2>, E3, not_a_Product>
{
  using type =
    add_t<multiply_t<E1, E3>, multiply_t<E2, E3>>;
};

template<Expression E1, Expression E2, Expression E3, Expression E4>
struct multiply<Sum<E1, E2>, Sum<E3, E4>, not_a_Product>
{
  using type =
    add_t<add_t<add_t<multiply_t<E1, E3>, multiply_t<E1, E4>>, multiply_t<E2, E3>>, multiply_t<E2, E4>>;
};

template<Expression E1, Expression E2>
struct make_power
{
  static_assert(is_symbol_v<E1> || is_power_v<E1>, "E1 must be a SymbolPowerType");
  static_assert(is_symbol_v<E2> || is_power_v<E2>, "E2 must be a SymbolPowerType");
  static_assert(E1::id_range.begin == E2::id_range.begin, "Can only combine powers of the same symbol!");

  using base = get_base_t<E1>;
  using new_exponent = typename add<get_exponent_t<E1>, get_exponent_t<E2>>::type;

  using type =
    /*if*/ std::conditional_t<is_constant_zero_v<new_exponent>,
      Constant<1, 1>,
    /*else if*/ std::conditional_t<is_constant_one_v<new_exponent>,
      base,
    /*else*/
      Power<base, new_exponent> >>;
};

template<int Enumerator1, int Denominator1, int Enumerator2, int Denominator2>
struct make_power<Constant<Enumerator1, Denominator1>, Constant<Enumerator2, Denominator2>>
{
  using type = multiply_t<Constant<Enumerator1, Denominator1>, Constant<Enumerator2, Denominator2>>;
};

} // namespace symbolic
