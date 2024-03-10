#pragma once

#include <type_traits>

namespace symbolic {

struct ExpressionTag {};

template<typename T>
concept Expression = std::is_base_of_v<ExpressionTag, T>;

// Helper class.
template<typename T>
struct DependentFalse : std::false_type
{
};

template<Expression E1, Expression E2>
struct expression_order_less
{
  consteval expression_order_less(E1 const&, E2 const&) { }

  consteval bool operator()() const
  {
    static_assert(DependentFalse<E1>::value, "Implement expression_order_less<>");
    return false;
  }
};

} // namespace symbolic
