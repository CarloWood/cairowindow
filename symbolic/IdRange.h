#pragma once

namespace symbolic {

template<int Begin, int End>
struct IdRange
{
  static constexpr int begin = Begin;
  static constexpr int end = End;
};

template<int begin1, int end1, int begin2, int end2>
consteval bool operator<(IdRange<begin1, end1>, IdRange<begin2, end2>)
{
  return end1 <= begin2;
}

} // namespace symbolic
