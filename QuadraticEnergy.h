#pragma once

#include "BezierCurve.h"
#include "symbolic/symbolic.h"

namespace cairowindow::autodiff {
using namespace symbolic;

class QuadraticEnergy : public BezierCurve
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
  static constexpr auto s_ = []() constexpr {
    return (v02_ + constant<2>() * z_ + a02_) ^ constant<1, 2>();
  }();
  static constexpr auto a03_ = []() constexpr {
    return a02_ * a0_;
  }();
  static constexpr auto za0pa03s_ = []() constexpr {
    return (z_ * a0_ + a03_) * s_;
  }();
  static constexpr auto za0v0_ = []() constexpr {
    return z_ * a0_ * v0_;
  }();
  static constexpr auto v02a02mz2_ = []() constexpr {
    return v02_ * a02_ - utils::square(z_);
  }();
  static constexpr auto zpa02pa0s_ = []() constexpr {
    return z_ + a02_ + a0_ * s_;
  }();
  static constexpr auto zpv0a0_ = []() constexpr {
    return z_ + v0_ * a0_;
  }();
  static constexpr auto the_log_ = []() constexpr {
    return log(zpa02pa0s_ / zpv0a0_);
  }();
  static constexpr auto the_enumerator_ = []() constexpr {
    return za0pa03s_ - za0v0_ + v02a02mz2_ * the_log_;
  }();
  static constexpr auto arc_length_ = []() constexpr {
    return the_enumerator_ / (constant<2>() * a03_);
  }();
#if 0
  static constexpr auto stretching_energy_ = []() constexpr {
    return utils::square(arc_length_);
  }();
#endif

 public:
  QuadraticEnergy(BezierCurve const orig) : BezierCurve(orig)
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

 public:
  double arc_length(double v0qa, double v1qa);
//  double stretching_energy(double v0qa, double v1qa);
};

} // namespace cairowindow::autodiff
