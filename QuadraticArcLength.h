#pragma once

#include "BezierCurve.h"
#include "symbolic/symbolic.h"
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
  static constexpr auto N1x_ = make_symbol();
  static constexpr auto N1y_ = make_symbol();
  static constexpr auto v0qa_ = make_symbol();
  static constexpr auto v1qa_ = make_symbol();

  static constexpr auto v0x_ = [](){
    auto v0_div_q1 = constant<2>() * sin(v1qa_) / sin(v1qa_ - v0qa_);
    auto t1 = cos(v0qa_) * Q1x_;
    auto t2 = sin(v0qa_) * N1x_;
    // multiply<
    //   Multiplication<
    //     Sin<Symbol<5>>,
    //     Exponentiation< Sin<Sum<Product<Constant<-1, 1>, Symbol<4>>, Symbol<5>>>, Constant<-1, 1>>
    //   >,
    //   Cos<Symbol<4>>,
    //   not_a_Product>
    auto s1 = v0_div_q1 * t1;
    auto s2 = v0_div_q1 * t2;
    auto r = s1 + s2;
    return r;
//    return v0_div_q1 * (cos(v0qa_) * Q1x_ + sin(v0qa_) * N1x_);
  }();

 public:
  QuadraticArcLength(BezierCurve const orig) : BezierCurve(orig)
  {
    Q1x_.register_name("Q1x");
    Q1y_.register_name("Q1y");
    N1x_.register_name("N1x");
    N1y_.register_name("N1y");
    v0qa_.register_name("alpha0");
    v1qa_.register_name("alpha1");
  }

  void test();

 public:
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
  double quadratic_arc_length();
};

} // namespace cairowindow::autodiff
