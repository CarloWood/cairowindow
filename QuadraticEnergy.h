#pragma once

#include "BezierCurve.h"
#include "symbolic/symbolic.h"

namespace cairowindow::autodiff {
using namespace symbolic;

class QuadraticEnergy : public BezierCurve
{
  Symbol const& Q1x_;
  Symbol const& Q1y_;
  Symbol const& N1x_;
  Symbol const& N1y_;

 public:
  Symbol const& v0qa_;
  Symbol const& v1qa_;

 private:
  Constant const& half = Constant::realize(1, 2);
  Constant const& two = Constant::realize(2);
  Constant const& three = Constant::realize(3);
  Constant const& four = Constant::realize(4);
  Constant const& five = Constant::realize(5);
  Constant const& six = Constant::realize(6);
  Constant const& c24 = Constant::realize(24);

  Function const& v0_div_q1_ = Function::realize("v0_div_q1", two * sin(v1qa_) / sin(v1qa_ - v0qa_));
  Function const& v0x_ = Function::realize("v0x", v0_div_q1_ * (cos(v0qa_) * Q1x_ + sin(v0qa_) * N1x_));
  Function const& v0y_ = Function::realize("v0y", v0_div_q1_ * (cos(v0qa_) * Q1y_ + sin(v0qa_) * N1y_));
  Function const& v02_ = Function::realize("v02", square(v0x_) + square(v0y_));
  Function const& A0x_ = Function::realize("A0x", two * (Q1x_ - v0x_));
  Function const& A0y_ = Function::realize("A0y", two * (Q1y_ - v0y_));
  Function const& a02_ = Function::realize("a02", square(A0x_) + square(A0y_));
  Function const& v0_ = Function::realize("v0", v02_ ^ half);
  Function const& a0_ = Function::realize("a0", a02_ ^ half);
  Function const& z_ = Function::realize("z", v0x_ * A0x_ + v0y_ * A0y_);
  Function const& s_ = Function::realize("s", (v02_ + two * z_ + a02_) ^ half);
  Function const& a03_ = Function::realize("a03", a02_ * a0_);
  Function const& za0pa03s_ = Function::realize("za0pa03s", (z_ * a0_ + a03_) * s_);
  Function const& za0v0_ = Function::realize("za0v0", z_ * a0_ * v0_);
  Function const& v02a02mz2_ = Function::realize("v02a02mz2", v02_ * a02_ - square(z_));
  Function const& zpa02pa0s_ = Function::realize("zpa02pa0s", z_ + a02_ + a0_ * s_);
  Function const& zpv0a0_ = Function::realize("zpv0a0", z_ + v0_ * a0_);
  Function const& the_log_ = Function::realize("the_log", log(zpa02pa0s_ / zpv0a0_));
  Function const& the_enumerator_ = Function::realize("the_enumerator", za0pa03s_ - za0v0_ + v02a02mz2_ * the_log_);
  Function const& arc_length_ = Function::realize("arc_length", the_enumerator_ / (two * a03_));
  Function const& stretching_energy_ = Function::realize("stretching_energy", square(arc_length_));

  Function const& dz_ = Function::realize("dz_", two * z_);
  Function const& cs_ = Function::realize("cs_", square(v0x_ * A0y_ - v0y_ * A0x_));
  Function const& a04_ = Function::realize("a04_", square(a02_));
  Function const& dz2_ = Function::realize("dz2_", square(dz_));
  Function const& e_ = Function::realize("e_", four * a02_ * v02_ - dz2_);
  Function const& e2_ = Function::realize("e2_", square(e_));
  Function const& se_ = Function::realize("se_", e_^half);
  Function const& i1_ = Function::realize("i1_",
      (-(dz_ + two * a02_) * (dz2_ - six * a02_ * dz_ - two * a02_ * (five * v02_ + three * a02_))) /
      square(a02_ + dz_ + v02_) + c24 * a04_ * atan((dz_ + two * a02_) / se_) / se_);
  Function const& i0_ = Function::realize("i0_", (-dz_ * (dz2_ - two * a02_ * five * v02_)) / square(v02_) + c24 * a04_ * atan(dz_ / se_) / se_);
  Function const& bending_energy_ = Function::realize("bending_energy_", cs_ * (i1_ - i0_) / (two * e2_));

 public:
  QuadraticEnergy(BezierCurve const orig, std::string name_prefix) : BezierCurve(orig),
      Q1x_(Symbol::realize(name_prefix + "Q1x")),
      Q1y_(Symbol::realize(name_prefix + "Q1y")),
      N1x_(Symbol::realize(name_prefix + "N1x")),
      N1y_(Symbol::realize(name_prefix + "N1y")),
      v0qa_(Symbol::realize(name_prefix + "v0qa")),
      v1qa_(Symbol::realize(name_prefix + "v1qa"))
  {
    // Q1 and N1 are constants (do not depend on v0qa or v1qa).
    Q1x_ = m_.coefficient[2].x() + m_.coefficient[1].x();
    Q1y_ = m_.coefficient[2].y() + m_.coefficient[1].y();
    N1x_ = -(m_.coefficient[2].y() + m_.coefficient[1].y());
    N1y_ = m_.coefficient[2].x() + m_.coefficient[1].x();
  }

 public:
  double arc_length(double v0qa, double v1qa);
  double stretching_energy(double v0qa, double v1qa);
  double bending_energy(double v0qa, double v1qa);

  Function const& arc_length() const { return arc_length_; }
  Function const& stretching_energy() const { return stretching_energy_; }
  Function const& bending_energy() const { return bending_energy_; }
};

} // namespace cairowindow::autodiff
