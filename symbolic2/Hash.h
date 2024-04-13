#pragma once

#include <boost/functional/hash.hpp>

namespace symbolic2 {

class Constant;
class Symbol;
class Sum;
class Product;
class Power;

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

template<>
struct Hash<Product>
{
  static constexpr uint64_t value = 0xA7D8AA34F62991AD;
};

template<>
struct Hash<Sum>
{
  static constexpr uint64_t value = 0x1DECA4B7CA07E34C;
};

template<>
struct Hash<Power>
{
  static constexpr uint64_t value = 0x1B156634A99A6A12;
};

} // namespace symbolic2
