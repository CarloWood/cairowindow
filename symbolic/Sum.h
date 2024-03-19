#pragma once

#include "IdRange.h"
#include "Expression.h"

namespace symbolic {

template<Expression E1, Expression E2>
class Sum : public ExpressionTag
{
 public:
  using arg1_type = E1;
  using arg2_type = E2;

  static constexpr precedence s_precedence = precedence::sum;
  static constexpr IdRange<E1::id_range.begin, E2::id_range.end> id_range{};

 private:
  E1 arg1_;
  E2 arg2_;

 public:
  inline constexpr Sum(E1 const& arg1, E2 const& arg2);

  E1 const& arg1() const { return arg1_; }
  E2 const& arg2() const { return arg2_; }

#ifdef SYMBOLIC_PRINTING
  void print_on(std::ostream& os) const
  {
    bool need_parens = needs_parens(arg1_.s_precedence, s_precedence, before);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " + ";
    need_parens = needs_parens(arg2_.s_precedence, s_precedence, after);
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<typename T>
struct is_sum : std::false_type { };

template<typename T>
constexpr bool is_sum_v = is_sum<T>::value;

template<Expression E1, Expression E2>
struct is_sum<Sum<E1, E2>> : std::true_type { };

template<typename E>
concept SumType = is_sum_v<E>;

template<Expression E1, Expression E2>
constexpr Sum<E1, E2>::Sum(E1 const& arg1, E2 const& arg2) : arg1_(arg1), arg2_(arg2)
{
  static_assert(!is_sum_v<E1>, "The first term of a Sum must not be a Sum itself.");
}

template<Expression E1, Expression E2, Expression E3, Expression E4>
constexpr bool is_same_expression(Sum<E1, E2> const& arg1, Sum<E3, E4> const& arg2)
{
  return std::is_same_v<E1, E3> && std::is_same_v<E2, E4>;
}

template<Expression E1, Expression E2>
auto operator+(E1 const& arg1, E2 const& arg2)
{
  if constexpr (!is_sum_v<E1> && !is_sum_v<E2>)
  {
    if constexpr (is_less_v<E1, E2>)
      return Sum{arg1, arg2};
    else if constexpr (is_less_v<E2, E1>)
      return Sum{arg2, arg1};
    else
      return constant<2>() * arg1;  //FIXME: add constants of arg1 and arg2.
  }
  else if constexpr (!is_sum_v<E1>)
  {
    if constexpr (is_less_v<E1, typename E2::arg1_type>)
      return Sum{arg1, arg2};
    else if constexpr (is_less_v<typename E2::arg1_type, E1>)
      return Sum{arg2.arg1(), arg1 + arg2.arg2()};
    else
      return constant<2>() * arg1 + arg2.arg2();  //FIXME: add constants of arg1 and arg2.
  }
  else if constexpr (!is_sum_v<E2>)
  {
    if constexpr (is_less_v<typename E1::arg1_type, E2>)
      return Sum{arg1.arg1(), arg2 + arg1.arg2()};
    else if constexpr (is_less_v<E2, typename E1::arg1_type>)
      return Sum{arg2, arg1};
    else
      return constant<2>() * arg1.arg1() + arg1.arg2();     //FIXME
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
      return constant<2>() * arg1.arg1() + arg1.arg2() + arg2.arg2();   // FIXME
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

} // namespace symbolic
