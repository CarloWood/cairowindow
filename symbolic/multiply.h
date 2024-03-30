#pragma once

#include "expression_traits.h"
#include "add.h"
#include "Product.h"

namespace symbolic {

template<Expression E1, Expression E2, IsProduct is_product>
struct multiply_unequals;

template<Expression E1, Expression E2>
struct multiply_unequals<E1, E2, isProduct>
{
  static consteval bool sanity_check()
  {
    if constexpr (is_product_v<E2>)
      return is_less_Product_v<E1, typename E2::arg1_type>;
    else
      return is_less_Product_v<E1, E2>;
  }

 private:
  static consteval auto eval()
  {
    static_assert(!is_product_v<E1>, "Only call multiply_unequals<..., isProduct> with a non-product arg1.");
    static_assert(sanity_check(), "Only call multiply_unequals for E1 < E2.");

    if constexpr (is_constant_zero_v<E1>)
      return Constant<0, 1>::instance();
    else if constexpr (is_constant_one_v<E1>)
      return E2::instance();
    else
      return Product<E1, E2>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2>
using make_product = typename multiply_unequals<E1, E2, isProduct>::type;

template<Expression E1, Expression E2>
struct make_power
{
 private:
  static consteval auto eval()
  {
    static_assert(is_symbol_v<E1> || is_power_v<E1>, "E1 must be a SymbolPowerType");
    static_assert(is_symbol_v<E2> || is_power_v<E2>, "E2 must be a SymbolPowerType");
    static_assert(E1::id_range.begin == E2::id_range.begin, "Can only combine powers of the same symbol!");

    using Base = get_base_t<E1>;
    using NewExponent = typename add<get_exponent_t<E1>, get_exponent_t<E2>>::type;

    if constexpr (is_constant_zero_v<NewExponent>)
      return Constant<1, 1>::instance();
    else if constexpr (is_constant_one_v<NewExponent>)
      return Base::instance();
    else
      return Power<Base, NewExponent>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2>
using make_power_t = typename make_power<E1, E2>::type;

template<int Enumerator1, int Denominator1, int Enumerator2, int Denominator2>
struct make_power<Constant<Enumerator1, Denominator1>, Constant<Enumerator2, Denominator2>>
{
  using type = multiply_t<Constant<Enumerator1, Denominator1>, Constant<Enumerator2, Denominator2>>;
};

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
      return make_product<E1, E2>::instance();
    else if constexpr (is_less_Product_v<E2, E1>)
      return make_product<E2, E1>::instance();
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
      return make_product<E1, Product<E2, E3>>::instance();
    else if constexpr (is_less_Product_v<E2, E1>)
      return make_product<E2, multiply_t<E1, E3>>::instance();
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
      return make_product<E1, multiply_t<E2, E3>>::instance();
    else if constexpr (is_less_Product_v<E3, E1>)
      return make_product<E3, Product<E1, E2>>::instance();
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
      return make_product<E1, multiply_t<E2, Product<E3, E4>>>::instance();
    else if constexpr (is_less_Product_v<E3, E1>)
      return make_product<E3, multiply_t<Product<E1, E2>, E4>>::instance();
    else
      return multiply_t<multiply_t<make_power_t<E1, E3>, E2>, E4>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2>
struct multiply_unequals<E1, E2, not_a_Product>
{
 private:
  template<typename T>
  static void ShowMeTheType()
  {
    static_assert(DependentFalse<T>::value, "WE GET HERE");
  }

  static consteval auto eval()
  {
    static_assert(!is_multiplication_v<E1>, "Only call multiply_unequals<..., not_a_Product> with a non-multiplication arg1.");
    static_assert(is_less_Multiplication_v<E1, E2>, "Only call multiply_unequals for E1 < E2.");

    // [with E1 = Constant<-1, 2>; E2 = Exponentiation<Sum<Symbol<10>, Symbol<11> >, Constant<3, 2> >]

    if constexpr (is_constant_zero_v<E1>)
      return Constant<0, 1>::instance();
    else if constexpr (is_constant_one_v<E1>)
      return E2::instance();
    else
    {
      using constant_factor1 = get_constant_factor_t<E1>;
      using constant_factor2 = get_constant_factor_t<E2>;
      using non_constant_factor2 = get_nonconstant_factor_t<E2>;
      using ConstantFactor = multiply_t<constant_factor1, constant_factor2>;

      static_assert(!is_constant_zero_v<ConstantFactor>, "Unexpected zero?!");

      if constexpr (is_constant_v<E1>)
      {
        if constexpr (is_constant_one_v<ConstantFactor>)
          return non_constant_factor2::instance();
        else
          return Multiplication<ConstantFactor, non_constant_factor2>::instance();
      }
      else
      {
        using non_constant_factor1 = get_nonconstant_factor_t<E1>;

        if constexpr (is_constant_one_v<ConstantFactor>)
          return Multiplication<non_constant_factor1, non_constant_factor2>::instance();
        else
        {
          using new_factor1 = multiply_t<ConstantFactor, non_constant_factor1>;
          if constexpr (is_product_v<new_factor1>)
            return Multiplication<new_factor1, non_constant_factor2>::instance();
          else
            return Multiplication<ConstantFactor, Multiplication<non_constant_factor1, non_constant_factor2>>::instance();
        }
      }
    }
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2>
using make_multiplication = typename multiply_unequals<E1, E2, not_a_Product>::type;

template<Expression E1, Expression E2>
struct make_exponentiation
{
 private:
  static consteval auto eval()
  {
#if 0
    using Base = get_base_t<E1>;
    using NewExponent = typename add<get_exponent_t<E1>, get_exponent_t<E2>>::type;

    if constexpr (is_constant_zero_v<NewExponent>)
      return Constant<1, 1>::instance();
    else if constexpr (is_constant_one_v<NewExponent>)
      return Base::instance();
    else
      return Power<Base, NewExponent>::instance();
#endif
    //FIXME
    return Multiplication<E1, E2>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2>
using make_exponentiation_t = typename make_exponentiation<E1, E2>::type;

template<Expression E1, Expression E2>
struct multiply<E1, E2, not_a_Product>
{
 private:
  static consteval auto eval()
  {
    if constexpr (is_less_Multiplication_v<E1, E2>)
      return make_multiplication<E1, E2>::instance();
    else if constexpr (is_less_Multiplication_v<E2, E1>)
      return make_multiplication<E2, E1>::instance();
    else
      return make_exponentiation_t<E1, E2>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2, Expression E3>
struct multiply<E1, Multiplication<E2, E3>, not_a_Product>
{
 private:
  static consteval auto eval()
  {
    if constexpr (is_less_Multiplication_v<E1, E2>)
      return make_multiplication<E1, Multiplication<E2, E3>>::instance();
    else if constexpr (is_less_Multiplication_v<E2, E1>)
      return make_multiplication<E2, multiply_t<E1, E3>>::instance();
    else
      return multiply_t<make_exponentiation_t<E1, E2>, E3>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2, Expression E3>
struct multiply<Multiplication<E1, E2>, E3, not_a_Product>
{
 private:
  static consteval auto eval()
  {
    if constexpr (is_less_Multiplication_v<E1, E3>)
      return make_multiplication<E1, multiply_t<E2, E3>>::instance();
    else if constexpr (is_less_Multiplication_v<E3, E1>)
      return make_multiplication<E3, Multiplication<E1, E2>>::instance();
    else
      return multiply_t<make_exponentiation_t<E3, E1>, E2>::instance();
  }

 public:
  using type = decltype(eval());
};

template<Expression E1, Expression E2, Expression E3, Expression E4>
struct multiply<Multiplication<E1, E2>, Multiplication<E3, E4>, not_a_Product>
{
  static constexpr bool E1_less_than_E3 = is_less_Multiplication_v<E1, E3>;
  static constexpr bool E3_less_than_E1 = is_less_Multiplication_v<E3, E1>;
  static constexpr bool E1_equal_E3 = !E1_less_than_E3 && !E3_less_than_E1;

 private:
  static consteval auto eval()
  {
    if constexpr (is_less_Multiplication_v<E1, E3>)
      return make_multiplication<E1, multiply_t<E2, Multiplication<E3, E4>>>::instance();
    else if constexpr (is_less_Multiplication_v<E3, E1>)
      return make_multiplication<E3, multiply_t<Multiplication<E1, E2>, E4>>::instance();
    else
      return multiply_t<multiply_t<make_exponentiation_t<E1, E3>, E2>, E4>::instance();
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

} // namespace symbolic
