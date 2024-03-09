#pragma once

#include "BezierCurve.h"
#include "utils/has_print_on.h"

namespace cairowindow::autodiff {
using utils::has_print_on::operator<<;

//----------------------------------------------------------------------------
// Classes for automatic differentiation.

struct AutoDiffTag { };
struct AutoDiffExpressionTag : public AutoDiffTag { };

template<typename T>
concept AutoDiff = std::is_base_of_v<AutoDiffTag, T>;

template<typename T>
concept AutoDiffExceptSymbol = std::is_base_of_v<AutoDiffExpressionTag, T>;

class Zero : public AutoDiffTag
{
 public:
  friend Zero operator-(Zero const&)
  {
    return {};
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << 0;
  }
#endif
};

template<typename T>
concept AutoDiffExceptZero = std::is_base_of_v<AutoDiffTag, T> && !std::is_same_v<T, Zero>;

class Constant : public AutoDiffTag
{
 private:
  double const value_;

 public:
  Constant(double value) : value_(value)
  {
    // Use class Zero.
    ASSERT(value != 0.0);
  }

  double value() const
  {
    return value_;
  }

  friend Constant operator-(Constant const& arg)
  {
    return {-arg.value()};
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << value_;
  }
#endif
};

template<AutoDiff T1>
class Negation;

template<uint32_t symbol_id_bit>
class Symbol : public AutoDiffTag
{
 private:
  static char const* s_name;
  static double s_value;

 public:
  Symbol(char const* name)
  {
    // Only instantiate any given Symbol once (did you use the same symbol_id_bit for two different symbols?).
    ASSERT(!s_name);
    s_name = name;
  }

  Symbol& operator=(double value)
  {
    s_value = value;
    return *this;
  }

  static double value()
  {
    return s_value;
  }

  static char const* name()
  {
    return s_name;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << s_name;
  }
#endif
};

template<uint32_t symbol_id_bit>
class DifferentiableSymbol : public Symbol<symbol_id_bit>
{
 public:
  DifferentiableSymbol(char const* name) : Symbol<symbol_id_bit>(name) { }

  using Symbol<symbol_id_bit>::operator=;
};

//static
template<uint32_t symbol_id_bit>
char const* Symbol<symbol_id_bit>::s_name;
//static
template<uint32_t symbol_id_bit>
double Symbol<symbol_id_bit>::s_value;

enum Symbols : uint32_t
{
  q1x,
  q1y,
  v0qa,
  v1qa
};

class QuadraticArcLength : public BezierCurve
{
 private:
  Symbol<q1x> Q1x{"Q1x"};
  Symbol<q1y> Q1y{"Q1y"};
  DifferentiableSymbol<v0qa> v0qa{"alpha0"};
  Symbol<v1qa> v1qa{"alpha1"};

 public:
  QuadraticArcLength(BezierCurve const orig) : BezierCurve(orig) { }

  void test();

 private:
  // Functions that returning intermediate values used to calculate quadratic_arc_length().
  double v0x(double v0qa, double v1qa);
  double v0y(double v0qa, double v1qa);
  double v02();
  double A0x();
  double A0y();
  double a02();
  double v0();
  double a0();
  double z();
  double s();
  double a03();
  double za0pa03s();
  double za0v0();
  double v02a02mz2();
  double zpa02pa0s();
  double zpv0a0();
  double the_log();
  double the_enumerator();
  double q_arc_length();
};

} // namespace cairowindow::autodiff
