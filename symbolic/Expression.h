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

// Define every is_less comparision as false, so that we
// only have to specialize the cases that are true.
template<Expression E1, Expression E2>
struct is_less : std::false_type { };

template<Expression E1, Expression E2>
constexpr bool is_less_v = is_less<E1, E2>::value;

} // namespace symbolic
