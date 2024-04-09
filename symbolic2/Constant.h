#pragma once

#include "Expression.h"
#include "Hash.h"
#include <numeric>
#include <ratio>
#include "debug.h"

namespace symbolic2 {

class Symbol;

class Constant : public Expression
{
 public:
  static Expression const& s_cached_zero;
  static Expression const& s_cached_one;

 private:
  int enumerator_;
  int denominator_;

  mutable uint64_t cached_hash_ = 0;

 private:
  Constant(int e, int d, int g) : enumerator_(e / g), denominator_(d / g) { }

  static int canonical_gcd(int e, int d)
  {
    int gcd = std::gcd(std::abs(e), std::abs(d));
    return d < 0 ? -gcd : gcd;
  }

  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Constant(int n) : enumerator_(n), denominator_(1) { }
  Constant(int e, int d) : Constant(e, d, canonical_gcd(e, d)) { }

  template<std::intmax_t Num, std::intmax_t Denom>
  Constant(std::ratio<Num, Denom>);

 public:
  static Constant const& realize(int n) { return static_cast<Constant const&>(get<Constant>(n)); }
  static Constant const& realize(int e, int d) { return static_cast<Constant const&>(get<Constant>(e, d)); }

  uint64_t hash() const override
  {
    if (cached_hash_ == 0)
    {
      cached_hash_ = hash_v<Constant>;
      boost::hash_combine(cached_hash_, enumerator_);
      boost::hash_combine(cached_hash_, denominator_);
    }

    return cached_hash_;
  }

  bool equals(Expression const& other) const override;

  double evaluate() const override
  {
    return static_cast<double>(enumerator_) / denominator_;
  }

  Expression const& differentiate(Symbol const&) const override
  {
    return s_cached_zero;
  }

#ifdef SYMBOLIC2_PRINTING
 public:
  void print_on(std::ostream& os) const override
  {
    os << enumerator_;
    if (denominator_ > 1)
      os << "/" << denominator_;
  }
#endif
};

template<std::intmax_t Num, std::intmax_t Denom>
Constant::Constant(std::ratio<Num, Denom>) : Constant(Num, Denom, canonical_gcd(Num, Denom))
{
}

} // namespace symbolic2
