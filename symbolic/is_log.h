#pragma once

#include "Expression.h"
#include <type_traits>

namespace symbolic {

template<Expression E>
class Log;

template<typename T>
struct is_log : std::false_type { };

template<typename T>
constexpr bool is_log_v = is_log<T>::value;

template<Expression E>
struct is_log<Log<E>> : std::true_type { };

template<typename T>
concept LogType = is_log_v<T>;

} // namespace symbolic
