#include "sys.h"
#include "QuadraticArcLength.h"
#include "symbolic/symbolic.h"

namespace cairowindow::autodiff {

// After calling BezierCurve::quadratic_from(double v0qa, double v1qa) (and that returned true),
// the following functions can be called to (redo) the calculations for the arc length;

// V₀_x.
double QuadraticArcLength::v0x(double v0qa, double v1qa)
{
#if 0
  // Q1 and N1 are constants (do not depend on v0qa or v1qa).
  double Q1x = m_.coefficient[2].x() + m_.coefficient[1].x();
  double N1x = -(m_.coefficient[2].y() + m_.coefficient[1].y());
  double v0_div_q1 = 2.0 * std::sin(v1qa) / std::sin(v1qa - v0qa);
  return v0_div_q1 * (std::cos(v0qa) * Q1x + std::sin(v0qa) * N1x);
#endif

  Q1x_ = m_.coefficient[2].x() + m_.coefficient[1].x();
  N1x_ = -(m_.coefficient[2].y() + m_.coefficient[1].y());
  v0qa_ = v0qa;
  v1qa_ = v1qa;

  return evaluate(v0x_);
}

// V₀_y.
double QuadraticArcLength::v0y(double v0qa, double v1qa)
{
  double Q1y = m_.coefficient[2].y() + m_.coefficient[1].y();
  double N1y = m_.coefficient[2].x() + m_.coefficient[1].x();
  double v0_div_q1 = 2.0 * std::sin(v1qa) / std::sin(v1qa - v0qa);
  return v0_div_q1 * (std::cos(v0qa) * Q1y + std::sin(v0qa) * N1y);
}

// |V₀|² = V₀_x² + V₀_y².
double QuadraticArcLength::v02()
{
  return V0().length_squared();
}

// A₀_x = 2 (Q1_x - V₀_x).
double QuadraticArcLength::A0x()
{
  double Q1x = m_.coefficient[2].x() + m_.coefficient[1].x();
  return 2.0 * (Q1x - V0().x());
}

// A₀_y = 2 (Q1_y - V₀_y).
double QuadraticArcLength::A0y()
{
  double Q1y = m_.coefficient[2].y() + m_.coefficient[1].y();
  return 2.0 * (Q1y - V0().y());
}

// |A₀|² = A₀_x² + A₀_y².
double QuadraticArcLength::a02()
{
  return A0().length_squared();
}

// v0 = |V₀| = √(V₀_x² + V₀_y²).
double QuadraticArcLength::v0()
{
  return std::sqrt(v02());
}

// a0 = |A₀| = √(A₀_x² + A₀_y²) = 2√((Q1_x - V₀_x)² + (Q1_y - V₀_y)²).
double QuadraticArcLength::a0()
{
  return std::sqrt(a02());
}

// z = V₀·A₀ = V₀_x A₀_x + V₀_y A₀_y = V₀_x (2 (Q1_x - V₀_x)) + V₀_y (2 (Q1_y - V₀_y)) =
//     2 (V₀_x Q1_x - V₀_x² + V₀_y Q1_y - V₀_y²).
double QuadraticArcLength::z()
{
  return V0().dot(A0());
}

double QuadraticArcLength::s()
{
  return std::sqrt(v02() + 2.0 * z() + a02());
}

// a03 = |A₀|³
double QuadraticArcLength::a03()
{
  return a02() * a0();
}

// (z * a0 + a03) * s
double QuadraticArcLength::za0pa03s()
{
  return (z() * a0() + a03()) * s();
}

double QuadraticArcLength::za0v0()
{
  return z() * a0() * v0();
}

double QuadraticArcLength::v02a02mz2()
{
  return v02() * a02() - utils::square(z());
}

double QuadraticArcLength::zpa02pa0s()
{
  return z() + a02() + a0() * s();
}

double QuadraticArcLength::zpv0a0()
{
  return z() + v0() * a0();
}

double QuadraticArcLength::the_log()
{
  return std::log(zpa02pa0s() / zpv0a0());
}

double QuadraticArcLength::the_enumerator()
{
  return za0pa03s() - za0v0() + v02a02mz2() * the_log();
}

double QuadraticArcLength::quadratic_arc_length()
{
  return the_enumerator() / (2.0 * a03());
}

// Test symbols.
namespace symbol {
static constexpr uint32_t x = 100;
static constexpr uint32_t y = 101;
} // namespace symbol

void QuadraticArcLength::test()
{
  auto zero1 = constant<0>();
  auto zero2 = constant<0>();

  Dout(dc::notice, "zero = " << zero1);
  Dout(dc::notice, "evaluate(zero) = " << evaluate(zero1));
  Dout(dc::notice, "-zero = " << (-zero1) << " = " << evaluate(-zero1));

  Dout(dc::notice, "zero + zero = " << (zero1 + zero2));
  Dout(dc::notice, "zero - zero = " << (zero1 - zero2));
  Dout(dc::notice, "zero * zero = " << (zero1 * zero2));

  auto constant1 = constant<42>();

  Dout(dc::notice, "constant = " << constant1);
  Dout(dc::notice, "evaluate(constant1) = " << evaluate(constant1));
  Dout(dc::notice, "-constant = " << (-constant1) << " = " << evaluate(-constant1));

  Dout(dc::notice, "zero + constant(" << evaluate(constant1) << ") = " << (zero1 + constant1));
  Dout(dc::notice, "zero - constant(" << evaluate(constant1) << ") = " << (zero1 - constant1));
  Dout(dc::notice, "zero * constant(" << evaluate(constant1) << ") = " << (zero1 * constant1));
  Dout(dc::notice, "zero / constant(" << evaluate(constant1) << ") = " << (zero1 / constant1));

  Dout(dc::notice, "constant(" << evaluate(constant1) << ") + zero = " << (constant1 + zero1));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") - zero = " << (constant1 - zero1));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") * zero = " << (constant1 * zero1));

  auto constant2 = constant<13>();

  Dout(dc::notice, "constant(" << evaluate(constant1) << ") + constant(" << evaluate(constant2) << ") = " << (constant1 + constant2));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") - constant(" << evaluate(constant2) << ") = " << (constant1 - constant2));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") * constant(" << evaluate(constant2) << ") = " << (constant1 * constant2));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") / constant(" << evaluate(constant2) << ") = " << (constant1 / constant2));

  auto x = make_symbol("x");
  x = 100;

  Dout(dc::notice, "x = " << x);
  Dout(dc::notice, "evaluate(x) = " << evaluate(x));

  Dout(dc::notice, "zero + x = " << (zero1 + x) << " = " << evaluate(zero1 + x));
  Dout(dc::notice, "zero - x = " << (zero1 - x) << " = " << evaluate(zero1 - x));
  Dout(dc::notice, "zero * x = " << (zero1 * x) << " = " << evaluate(zero1 * x));
  Dout(dc::notice, "zero / x = " << (zero1 / x) << " = " << evaluate(zero1 / x));

  Dout(dc::notice, "x + zero = " << (x + zero1) << " = " << evaluate(x + zero1));
  Dout(dc::notice, "x - zero = " << (x - zero1) << " = " << evaluate(x - zero1));
  Dout(dc::notice, "x * zero = " << (x * zero1) << " = " << evaluate(x * zero1));

  Dout(dc::notice, "constant(" << evaluate(constant1) << ") + x = " << (constant1 + x) << " = " << evaluate(constant1 + x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") - x = " << (constant1 - x) << " = " << evaluate(constant1 - x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") * x = " << (constant1 * x) << " = " << evaluate(constant1 * x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") / x = " << (constant1 / x) << " = " << evaluate(constant1 / x));

  Dout(dc::notice, "x + constant(" << evaluate(constant1) << ") = " << (x + constant1) << " = " << evaluate(x + constant1));
  Dout(dc::notice, "x - constant(" << evaluate(constant1) << ") = " << (x - constant1) << " = " << evaluate(x - constant1));
  Dout(dc::notice, "x * constant(" << evaluate(constant1) << ") = " << (x * constant1) << " = " << evaluate(x * constant1));
  Dout(dc::notice, "x / constant(" << evaluate(constant1) << ") = " << (x / constant1) << " = " << evaluate(x / constant1));

  Dout(dc::notice, "-x = " << (-x) << " = " << evaluate(-x));

  Dout(dc::notice, "zero + -x = " << (zero1 + -x) << " = " << evaluate(zero1 + -x));
  Dout(dc::notice, "zero - -x = " << (zero1 - -x) << " = " << evaluate(zero1 - -x));
  Dout(dc::notice, "zero * -x = " << (zero1 * -x) << " = " << evaluate(zero1 * -x));
  Dout(dc::notice, "zero / -x = " << (zero1 / -x) << " = " << evaluate(zero1 / -x));

  Dout(dc::notice, "-x + zero = " << (-x + zero1) << " = " << evaluate(-x + zero1));
  Dout(dc::notice, "-x - zero = " << (-x - zero1) << " = " << evaluate(-x - zero1));
  Dout(dc::notice, "-x * zero = " << (-x * zero1) << " = " << evaluate(-x * zero1));

  Dout(dc::notice, "constant(" << evaluate(constant1) << ") + -x = " << (constant1 + -x) << " = " << evaluate(constant1 + -x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") - -x = " << (constant1 - -x) << " = " << evaluate(constant1 - -x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") * -x = " << (constant1 * -x) << " = " << evaluate(constant1 * -x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") / -x = " << (constant1 / -x) << " = " << evaluate(constant1 / -x));

  Dout(dc::notice, "-x + constant(" << evaluate(constant1) << ") = " << (-x + constant1) << " = " << evaluate(-x + constant1));
  Dout(dc::notice, "-x - constant(" << evaluate(constant1) << ") = " << (-x - constant1) << " = " << evaluate(-x - constant1));
  Dout(dc::notice, "-x * constant(" << evaluate(constant1) << ") = " << (-x * constant1) << " = " << evaluate(-x * constant1));
  Dout(dc::notice, "-x / constant(" << evaluate(constant1) << ") = " << (-x / constant1) << " = " << evaluate(-x / constant1));

  auto y = make_symbol("y");
  y = 200;

  Dout(dc::notice, "x + y = " << (x + y) << " = " << evaluate(x + y));
  Dout(dc::notice, "x - y = " << (x - y) << " = " << evaluate(x - y));
  Dout(dc::notice, "x * y = " << (x * y) << " = " << evaluate(x * y));
  Dout(dc::notice, "x / y = " << (x / y) << " = " << evaluate(x / y));

  Dout(dc::notice, "-x + y = " << (-x + y) << " = " << evaluate(-x + y));
  Dout(dc::notice, "-x - y = " << (-x - y) << " = " << evaluate(-x - y));
  Dout(dc::notice, "-x * y = " << (-x * y) << " = " << evaluate(-x * y));
  Dout(dc::notice, "-x / y = " << (-x / y) << " = " << evaluate(-x / y));

  Dout(dc::notice, "x + -y = " << (x + -y) << " = " << evaluate(x + -y));
  Dout(dc::notice, "x - -y = " << (x - -y) << " = " << evaluate(x - -y));
  Dout(dc::notice, "x * -y = " << (x * -y) << " = " << evaluate(x * -y));
  Dout(dc::notice, "x / -y = " << (x / -y) << " = " << evaluate(x / -y));

  Dout(dc::notice, "-x + -y = " << (-x + -y) << " = " << evaluate(-x + -y));
  Dout(dc::notice, "-x - -y = " << (-x - -y) << " = " << evaluate(-x - -y));
  Dout(dc::notice, "-x * -y = " << (-x * -y) << " = " << evaluate(-x * -y));
  Dout(dc::notice, "-x / -y = " << (-x / -y) << " = " << evaluate(-x / -y));

  auto formula1 = (x + y)^constant<3, 2>();
  Dout(dc::notice, "(x + y)^(3/2) = " << formula1 << " = " << evaluate(formula1));
  auto formula2 = formula1 * (constant<7>() * x - constant<1, 7>() * y);
  Dout(dc::notice, "(x + y)^(3/2) * (7 * x - 1/7 * y) = " << formula2 << " = " << evaluate(formula2));

  Dout(dc::notice, differentiate(constant1, x));

  Dout(dc::notice, "∂x/∂x = " << differentiate(x, x));
  Dout(dc::notice, "∂y/∂x = " << differentiate(y, x));

  Dout(dc::notice, "∂-x/∂x = " << differentiate(-x, x));
  Dout(dc::notice, "∂(x + y + 42)/∂x = " << differentiate(x + y + constant1, x));
  Dout(dc::notice, "∂(x + y - 42)/∂x = " << differentiate(x + y - constant1, x));
  Dout(dc::notice, "∂(y - x)/∂x = " << differentiate(y - x, x));

  Dout(dc::notice, "∂(x * x * x)/∂x = " << differentiate(x * x * x, x));

  Dout(dc::notice, "∂((y - x^2)^5)/∂x = " << differentiate((y - (x^constant<2>()))^constant<5>(), x));
}

} // namespace cairowindow::autodiff
