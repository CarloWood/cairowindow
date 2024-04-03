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

  static constexpr auto v0_div_q1_ = []() constexpr {
    return constant<2>() * sin(v1qa_) / sin(v1qa_ - v0qa_);
  }();
  static constexpr auto v0x_ = []() constexpr {
    return v0_div_q1_ * (cos(v0qa_) * Q1x_ + sin(v0qa_) * N1x_);
  }();
  static constexpr auto v0y_ = []() constexpr {
    return v0_div_q1_ * (cos(v0qa_) * Q1y_ + sin(v0qa_) * N1y_);
  }();
  static constexpr auto v02_ = []() constexpr {
    return utils::square(v0x_) + utils::square(v0y_);
  }();
  static constexpr auto A0x_ = []() constexpr {
    return constant<2>() * (Q1x_ - v0x_);
  }();
  static constexpr auto A0y_ = []() constexpr {
    return constant<2>() * (Q1y_ - v0y_);
  }();
  static constexpr auto a02_ = []() constexpr {
    return utils::square(A0x_) + utils::square(A0y_);
  }();
  static constexpr auto v0_ = []() constexpr {
    return v02_ ^ constant<1, 2>();
  }();
  static constexpr auto a0_ = []() constexpr {
    return a02_ ^ constant<1, 2>();
  }();
  static constexpr auto z_ = []() constexpr {
    return v0x_ * A0x_ + v0y_ * A0y_;
  }();
 public:
  static constexpr auto s_ = []() constexpr {
    return (v02_ + constant<2>() * z_ + a02_) ^ constant<1, 2>();
  }();
  static constexpr auto a03_ = []() constexpr {
    return a02_ * a0_;
  }();
  static constexpr auto za0pa03s_ = []() constexpr {
    return (z_ * a0_ + a03_) * s_;
  }();

 public:
  QuadraticArcLength(BezierCurve const orig) : BezierCurve(orig)
  {
    Q1x_.register_name("Q1x");
    Q1y_.register_name("Q1y");
    N1x_.register_name("N1x");
    N1y_.register_name("N1y");

    // Q1 and N1 are constants (do not depend on v0qa or v1qa).
    Q1x_ = m_.coefficient[2].x() + m_.coefficient[1].x();
    Q1y_ = m_.coefficient[2].y() + m_.coefficient[1].y();
    N1x_ = -(m_.coefficient[2].y() + m_.coefficient[1].y());
    N1y_ = m_.coefficient[2].x() + m_.coefficient[1].x();

    v0qa_.register_name("alpha0");
    v1qa_.register_name("alpha1");
  }

  void test();

 public:
  double test_v0qa_;
  double test_v1qa_;

  // Functions that returning intermediate values used to calculate quadratic_arc_length().
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
  double quadratic_arc_length(double v0qa, double v1qa);
};

} // namespace cairowindow::autodiff
