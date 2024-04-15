#pragma once

#include <boost/functional/hash.hpp>

namespace symbolic {

class Constant;
class Symbol;
class Sum;
class Product;
class Power;
class Sin;
class Cos;
class Atan;
class Log;
class Function;

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

template<>
struct Hash<Sin>
{
  static constexpr uint64_t value = 0x8B92C1E8AB8BF76B;
};

template<>
struct Hash<Cos>
{
  static constexpr uint64_t value = 0x514EC2A8F6B3F521;
};

template<>
struct Hash<Atan>
{
  static constexpr uint64_t value = 0xBD0BF25373FC829F;
};

template<>
struct Hash<Log>
{
  static constexpr uint64_t value = 0x3F3EBCCAA947563D;
};

template<>
struct Hash<Function>
{
  static constexpr uint64_t value = 0x0904BFE81CCFDE47;
};

} // namespace symbolic
