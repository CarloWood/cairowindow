#pragma once

#include "IdRange.h"
#include "Expression.h"
#include <algorithm>

namespace symbolic {

template<Expression E1, Expression E2>
class Sum : public ExpressionTag
{
 public:
  using arg1_type = E1;
  using arg2_type = E2;

  static constexpr precedence s_precedence = precedence::sum;
  static constexpr IdRange<std::min(E1::id_range.begin, E2::id_range.begin), std::max(E1::id_range.end, E2::id_range.end)> id_range{};

 private:
  E1 arg1_;
  E2 arg2_;

 public:
  inline constexpr Sum(E1 const& arg1, E2 const& arg2);

  E1 const& arg1() const { return arg1_; }
  E2 const& arg2() const { return arg2_; }

#ifdef SYMBOLIC_PRINTING
  void print_on(std::ostream& os) const;
#endif
};

template<typename T>
struct is_sum : std::false_type { };

template<typename T>
constexpr bool is_sum_v = is_sum<T>::value;

template<Expression E1, Expression E2>
struct is_sum<Sum<E1, E2>> : std::true_type { };

#ifdef SYMBOLIC_PRINTING
template<Expression E1, Expression E2>
void Sum<E1, E2>::print_on(std::ostream& os) const
{
  bool need_parens = needs_parens(arg1_.s_precedence, s_precedence, before);
  if (need_parens)
    os << '(';
  arg1_.print_on(os);
  if (need_parens)
    os << ')';
  if constexpr (is_product_v<E2>)
  {
    if constexpr (is_constant_v<typename E2::arg1_type>)
    {
      if constexpr (E2::arg1_type::is_less_than_zero())
      {
        os << " - ";
        auto term = -arg2_;
        need_parens = needs_parens(term.s_precedence, precedence::difference, after);
        if (need_parens)
          os << '(';
        term.print_on(os);
        if (need_parens)
          os << ')';
        return;
      }
    }
  }
  else if constexpr (is_sum_v<E2>)
  {
    if constexpr (is_product_v<typename E2::arg1_type>)
    {
      if constexpr (is_constant_v<typename E2::arg1_type::arg1_type>)
      {
        if constexpr (E2::arg1_type::arg1_type::is_less_than_zero())
        {
          os << " - ";
          (-arg2_.arg1() + arg2_.arg2()).print_on(os);
          return;
        }
      }
    }
  }
  os << " + ";
  need_parens = needs_parens(arg2_.s_precedence, s_precedence, after);
  if (need_parens)
    os << '(';
  arg2_.print_on(os);
  if (need_parens)
    os << ')';
}
#endif

template<typename E>
concept SumType = is_sum_v<E>;

template<Expression E1, Expression E2>
constexpr Sum<E1, E2>::Sum(E1 const& arg1, E2 const& arg2) : arg1_(arg1), arg2_(arg2)
{
  static_assert(!is_sum_v<E1>, "The first term of a Sum must not be a Sum itself.");
  static_assert(!is_constant_v<E2>, "The second term of a Sum is not allowed to be a constant.");
  if constexpr (is_constant_v<E1>)
    static_assert(!E1::is_zero(), "Don't construct a Sum with a zero constant.");
}

template<Expression E1, Expression E2, Expression E3, Expression E4>
constexpr bool is_same_expression(Sum<E1, E2> const& arg1, Sum<E3, E4> const& arg2)
{
  return std::is_same_v<E1, E3> && std::is_same_v<E2, E4>;
}

template<Expression E>
struct get_constant_factor
{
  using type = Constant<1, 1>;
};

template<Expression E>
using get_constant_factor_t = typename get_constant_factor<E>::type;

template<int Enumerator1, int Denominator1, Expression E2>
struct get_constant_factor<Product<Constant<Enumerator1, Denominator1>, E2>>
{
  using type = Constant<Enumerator1, Denominator1>;
};

template<Expression E>
auto const& get_nonconstant_factor(E const& arg)
{
  return arg;
}

template<int Enumerator1, int Denominator1, Expression E2>
auto const& get_nonconstant_factor(Product<Constant<Enumerator1, Denominator1>, E2> const& arg)
{
  return arg.arg2();
}

template<Expression E1, Expression E2>
constexpr auto add_equals(E1 const& arg1, E2 const& arg2)
{
  static_assert(!is_sum_v<E1>, "The first term of add_equals must not be a Sum.");

  if constexpr (is_constant_v<E1>)
  {
    static_assert(is_constant_v<E2>, "!");
    return constant<E1::s_enumerator * E2::s_denominator + E2::s_enumerator * E1::s_denominator, E1::s_denominator * E2::s_denominator>();
  }
  else
  {
    using constant_factor1 = get_constant_factor_t<E1>;
    using constant_factor2 = get_constant_factor_t<E2>;
    auto const& nonconstant_factor = get_nonconstant_factor(arg1);
    static_assert(std::is_same_v<decltype(nonconstant_factor), decltype(get_nonconstant_factor(arg2))>, "Expected the same non-constant type.");

    constexpr auto constant_factor = constant<constant_factor1::s_enumerator * constant_factor2::s_denominator +
      constant_factor2::s_enumerator * constant_factor1::s_denominator, constant_factor1::s_denominator * constant_factor2::s_denominator>();

    if constexpr (constant_factor.is_zero())
      return constant<0>();
    else if constexpr (constant_factor.is_one())
      return nonconstant_factor;
    else
      return Product{constant_factor, nonconstant_factor};
  }
}

template<Expression E1, Expression E2>
struct is_less : std::false_type { };

template<Expression E1, Expression E2>
constexpr bool is_less_v = is_less<E1, E2>::value;

// A Constant is (only) less than a non-constant.
template<int Enumerator1, int Denominator1, Expression E2>
requires (!is_constant_v<E2>)
struct is_less<Constant<Enumerator1, Denominator1>, E2> : std::true_type { };

// Compare two symbols.
template<int Id1, int Id2>
requires (Id1 < Id2)
struct is_less<Symbol<Id1>, Symbol<Id2>> : std::true_type { };

// Compare a symbol with something not a constant or symbol.
template<int Id1, Expression E2>
requires (!is_constant_v<E2> && !is_symbol_v<E2> && !is_product_v<E2>)
struct is_less<Symbol<Id1>, E2> : std::true_type { };

template<int Id1, Expression E2, Expression E3>
requires (!is_constant_v<E2> || is_less_v<Symbol<Id1>, E3>)
struct is_less<Symbol<Id1>, Product<E2, E3>> : std::true_type { };

template<Expression E1, Expression E2, int Id3>
requires (is_constant_v<E1> && is_less_v<E2, Symbol<Id3>>)
struct is_less<Product<E1, E2>, Symbol<Id3>> : std::true_type { };

template<SymbolType E1, int Enumerator1, int Denominator1, SymbolType E2, int Enumerator2, int Denominator2>
requires (is_less_v<E1, E2> || (!is_less_v<E2, E1> && Enumerator1 * Denominator2 < Enumerator2 * Denominator1))
struct is_less<Power<E1, Enumerator1, Denominator1>, Power<E2, Enumerator2, Denominator2>> : std::true_type { };

template<SymbolType E1, int Enumerator1, int Denominator1, Expression E2>
requires (!is_constant_v<E2> && !is_symbol_v<E2> && !is_power_v<E2> && !is_product_v<E2>)
struct is_less<Power<E1, Enumerator1, Denominator1>, E2> : std::true_type { };

template<SymbolType E1, int Enumerator1, int Denominator1, Expression E2, Expression E3>
requires (!is_constant_v<E2> || is_less_v<Power<E1, Enumerator1, Denominator1>, E3>)
struct is_less<Power<E1, Enumerator1, Denominator1>, Product<E2, E3>> : std::true_type { };

template<Expression E1, Expression E2, SymbolType S, int Enumerator, int Denominator>
requires (is_constant_v<E1> && is_less_v<E2, Power<S, Enumerator, Denominator>>)
struct is_less<Product<E1, E2>, Power<S, Enumerator, Denominator>> : std::true_type { };

template<Expression E1, Expression E2, Expression E3, Expression E4>
requires (
    (!is_constant_v<E1> && !is_constant_v<E3> &&
     (is_less_v<E1, E3> || (!is_less_v<E3, E1> && is_less_v<E2, E4>))) ||
    (!is_constant_v<E1> && is_constant_v<E3> && is_less_v<Product<E1, E2>, E4>) ||
    (is_constant_v<E1> && !is_constant_v<E3> && is_less_v<E2, Product<E3, E4>>) ||
    (is_constant_v<E1> && is_constant_v<E3> && is_less_v<E2, E4>))
struct is_less<Product<E1, E2>, Product<E3, E4>> : std::true_type { };

template<Expression E1, Expression E2, Expression E3>
requires (!is_constant_v<E3> && !is_symbol_v<E3> && !is_power_v<E3> && !is_product_v<E3>)
struct is_less<Product<E1, E2>, E3> : std::true_type { };

template<Expression E1, Expression E2>
auto operator+(E1 const& arg1, E2 const& arg2)
{
  if constexpr (!is_sum_v<E1> && !is_sum_v<E2>)
  {
    if constexpr (is_less_v<E1, E2>)
    {
      if constexpr (is_constant_v<E1>)
      {
        if constexpr (E1::is_zero())
          return arg2;
        else
          return Sum{arg1, arg2};
      }
      else
        return Sum{arg1, arg2};
    }
    else if constexpr (is_less_v<E2, E1>)
    {
      if constexpr (is_constant_v<E2>)
      {
        if constexpr (E2::is_zero())
          return arg1;
        else
          return Sum{arg2, arg1};
      }
      else
        return Sum{arg2, arg1};
    }
    else
      return add_equals(arg1, arg2);
  }
  else if constexpr (!is_sum_v<E1>)
  {
    if constexpr (is_less_v<E1, typename E2::arg1_type>)
    {
      if constexpr (is_constant_v<E1>)
      {
        if constexpr (E1::is_zero())
          return arg2;
        else
          return Sum{arg1, arg2};
      }
      else
        return Sum{arg1, arg2};
    }
    else if constexpr (is_less_v<typename E2::arg1_type, E1>)
      return Sum{arg2.arg1(), arg1 + arg2.arg2()};
    else
      return add_equals(arg1, arg2.arg1()) + arg2.arg2();
  }
  else if constexpr (!is_sum_v<E2>)
  {
    if constexpr (is_less_v<typename E1::arg1_type, E2>)
      return Sum{arg1.arg1(), arg2 + arg1.arg2()};
    else if constexpr (is_less_v<E2, typename E1::arg1_type>)
      return Sum{arg2, arg1};
    else
      return add_equals(arg1.arg1(), arg2) + arg1.arg2();
  }
  else
  {
    if constexpr (is_less_v<typename E1::arg1_type, typename E2::arg1_type>)
    {
      return Sum{arg1.arg1(), arg1.arg2() + arg2};
    }
    else if constexpr (is_less_v<typename E2::arg1_type, typename E1::arg1_type>)
    {
      return Sum{arg2.arg1(), arg1 + arg2.arg2()};
    }
    else
      return add_equals(arg1.arg1(), arg2.arg1()) + arg1.arg2() + arg2.arg2();
  }
}

template<Expression E1, Expression E2>
constexpr auto inverse(Sum<E1, E2> sum)
{
  return Power<Sum<E1, E2>, -1, 1>{sum};
}

template<Expression E1, Expression E2, Expression E3, Expression E4>
requires (is_less_v<E1, E3> || (!is_less_v<E3, E1> && is_less_v<E2, E4>))
struct is_less<Sum<E1, E2>, Sum<E3, E4>> : std::true_type { };

template<Expression E1, Expression E2, Expression E3>
requires (!is_constant_v<E3> && !is_symbol_v<E3> && !is_power_v<E3> && !is_product_v<E3>)
struct is_less<Sum<E1, E2>, E3> : std::true_type { };

template<Expression E1, Expression E2>
auto operator-(E1 const& arg1, E2 const& arg2)
{
  return arg1 + -arg2;
}

template<Expression E1, Expression E2, Expression E3>
requires (!is_sum_v<E1>)
auto operator*(E1 const& arg1, Sum<E2, E3> const& arg2)
{
  return arg1 * arg2.arg1() + arg1 * arg2.arg2();
}

template<Expression E1, Expression E2, Expression E3>
requires (!is_sum_v<E3>)
auto operator*(Sum<E1, E2> const& arg1, E3 const& arg2)
{
  return arg1.arg1() * arg2 + arg1.arg2() * arg2;
}

template<Expression E1, Expression E2, Expression E3, Expression E4>
auto operator*(Sum<E1, E2> const& arg1, Sum<E3, E4> const& arg2)
{
  return arg1.arg1() * arg2.arg1() + arg1.arg1() * arg2.arg2() + arg1.arg2() * arg2.arg1() + arg1.arg2() * arg2.arg2();
}

} // namespace symbolic
