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

} // namespace symbolic
