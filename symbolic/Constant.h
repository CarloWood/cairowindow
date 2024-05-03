#pragma once

#include "Expression.h"
#include "Hash.h"
#include <boost/functional/hash.hpp>
#include <numeric>
#include <ratio>
#include "debug.h"

namespace symbolic {

class Symbol;

class Constant : public Expression
{
 public:
  static Constant const& s_cached_zero;
  static Constant const& s_cached_one;
  static Constant const& s_cached_minus_one;

  int const enumerator_;
  int const denominator_;

 private:
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

  Precedence precedence() const override final { return denominator_ != 1 ? Precedence::ratio : enumerator_ < 0 ? Precedence::negation : Precedence::constant; }

  Expression const& get_nonconstant_factor() const override final { return s_cached_one; }

 public:
  static Constant const& realize(int n)
  {
    DoutEntering(dc::symbolic, "Constant::realize(" << n << ")");
    return static_cast<Constant const&>(get<Constant>(n));
  }
  static Constant const& realize(int e, int d)
  {
    DoutEntering(dc::symbolic, "Constant::realize(" << e << ", " << d << ")");
    return static_cast<Constant const&>(get<Constant>(e, d));
  }

  static bool is_zero(Expression const& arg) { return &arg == &s_cached_zero; }
  static bool is_one(Expression const& arg) { return &arg == &s_cached_one; }
  static bool is_minus_one(Expression const& arg) { return &arg == &s_cached_minus_one; }

  bool is_negative() const { return enumerator_ < 0; }

  ExpressionType type() const override final { return constantT; }

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

  Expression const& derivative(Symbol const&) const override
  {
    return s_cached_zero;
  }

  bool operator<(Constant const& other) const { return enumerator_ * other.denominator_ < other.enumerator_ * denominator_; }

  friend Constant const& operator+(Constant const& arg1, Constant const& arg2)
  {
    return realize(arg1.enumerator_ * arg2.denominator_ + arg2.enumerator_ * arg1.denominator_, arg1.denominator_ * arg2.denominator_);
  }

  friend Constant const& operator-(Constant const& arg1, Constant const& arg2)
  {
    return realize(arg1.enumerator_ * arg2.denominator_ - arg2.enumerator_ * arg1.denominator_, arg1.denominator_ * arg2.denominator_);
  }

  friend Constant const& operator*(Constant const& arg1, Constant const& arg2)
  {
    return realize(arg1.enumerator_ * arg2.enumerator_, arg1.denominator_ * arg2.denominator_);
  }

  friend Constant const& operator/(Constant const& arg1, Constant const& arg2)
  {
    return realize(arg1.enumerator_ * arg2.denominator_, arg1.denominator_ * arg2.enumerator_);
  }

#ifdef SYMBOLIC_PRINTING
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

} // namespace symbolic
