#pragma once

#include "utils/type_id_hash.h"

namespace symbolic {

template<typename T>
constexpr uint64_t hash_v = utils::type_id_hash<T>();

} // namespace symbolic
