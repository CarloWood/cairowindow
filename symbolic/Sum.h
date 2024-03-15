#pragma once

#include "IdRange.h"
#include "Expression.h"

namespace symbolic {

template<Expression E1, Expression E2>
class Sum : public ExpressionTag
{
 public:
  static constexpr precedence s_precedence = precedence::sum;
  static constexpr IdRange<E1::id_range.begin, E2::id_range.end> id_range{};

 private:
  E1 arg1_;
  E2 arg2_;

 public:
  inline constexpr Sum(E1 const& arg1, E2 const& arg2);

  constexpr E1 const& arg1() const { return arg1_; }
  constexpr E2 const& arg2() const { return arg2_; }

#ifdef SYMBOLIC_PRINTING
  constexpr bool needs_parens(before_or_after position, precedence prec) const { return prec <= s_precedence; }
  constexpr bool needs_parens(precedence prec) const { return prec <= s_precedence; }

  void print_on(std::ostream& os) const
  {
    bool need_parens = arg1_.needs_parens(before, s_precedence);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " + ";
    need_parens = arg2_.needs_parens(after, s_precedence);
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
  return Sum{arg1, arg2};
}

template<Expression E1, Expression E2>
constexpr auto inverse(Sum<E1, E2> sum)
{
  return Power<Sum<E1, E2>, -1, 1>{sum};
}

} // namespace symbolic
