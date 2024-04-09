#pragma once

#include <boost/functional/hash.hpp>

namespace symbolic2 {

class Constant;
class Symbol;

template<typename T>
struct Hash;

template<typename T>
constexpr uint64_t hash_v = Hash<T>::value;

template<>
struct Hash<Constant>
{
  static constexpr uint64_t value = 0x5C2D67C9A0D3CF0C;
};

template<>
struct Hash<Symbol>
{
  static constexpr uint64_t value = 0xB458550B80FF25DA;
};

} // namespace symbolic2
