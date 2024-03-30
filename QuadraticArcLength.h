#pragma once

#include "BezierCurve.h"
#include "symbolic/Symbol.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace cairowindow::autodiff {
using namespace symbolic;
#ifdef CWDEBUG
//using utils::has_print_on::operator<<;
#endif

class QuadraticArcLength : public BezierCurve
{
 private:
  static constexpr auto Q1x_ = make_symbol();
  static constexpr auto Q1y_ = make_symbol();
  static constexpr auto v0qa_ = make_symbol();
  static constexpr auto v1qa_ = make_symbol();

 public:
  QuadraticArcLength(BezierCurve const orig) : BezierCurve(orig)
  {
    Q1x_.register_name("Q1x");
    Q1y_.register_name("Q1y");
    v0qa_.register_name("alpha0");
    v1qa_.register_name("alpha1");
  }

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
