#pragma once

#include <numeric>

namespace symbolic {

template<int Enumerator, int Denominator>
class Constant;

template<int Enumerator, int Denominator>
struct canonical_constant
{
  static constexpr int abs_enumerator = Enumerator < 0 ? -Enumerator : Enumerator;
  static constexpr int abs_denominator = Denominator < 0 ? -Denominator : Denominator;
  static constexpr int signed_enumerator = Denominator < 0 ? -Enumerator : Enumerator;
  static constexpr int gcd = std::gcd(abs_enumerator, abs_denominator);
  using type = Constant<signed_enumerator / gcd, abs_denominator / gcd>;
};

template<int Enumerator, int Denominator>
using canonical_constant_t = typename canonical_constant<Enumerator, Denominator>::type;

} // namespace symbolic
