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
  // If the symbol id ranges do not overlap, rearrange them if necessary
  // so that the symbols with the smallest id are on the left.
  if constexpr (E2::id_range < E1::id_range)
  {
    // However, first term is not allowed to be a Sum itself.
    // Therefore, if the second term is a Sum then we can't just swap the two arguments.
    if constexpr (is_sum_v<E2>)
      return Sum{arg2.arg1(), arg2.arg2() + arg1};      // (k + L) + (a + B) --> a + (B + (k + L)).
                                                        // Note that the range of B is guaranteed to be less than that of k + L.
    else
      return Sum{arg2, arg1};                           // arg2 is not a Sum, so simply swapping is allowed.
  }
  else if constexpr (is_sum_v<E1>)
    return arg1.arg1() + (arg1.arg2() + arg2);          // (a + b) + c = a + (b + c).
  else
    return Sum{arg1, arg2};
}

template<Expression E1, Expression E2>
constexpr auto inverse(Sum<E1, E2> sum)
{
  return Power<Sum<E1, E2>, -1, 1>{sum};
}

} // namespace symbolic
