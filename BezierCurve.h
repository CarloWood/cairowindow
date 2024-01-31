#pragma once

#include "Vector.h"
#include "Rectangle.h"
#include "utils/has_print_on.h"
#include <array>
#ifdef CWDEBUG
#include "debug.h"
#endif

namespace cairowindow {
using utils::has_print_on::operator<<;

struct BezierCurveMatrix
{
  std::array<Vector, 4> coefficient;
};

// The cubic Bezier curve is the 2D parametric function:
//
//   P(t) = (1-t)³ P₀ + 3(1-t)²t C₁ + 3(1-t)t² C₂ + t³ P₁
//
//        = (1 - 3t + 3t² - t³) P₀ + (3t - 6t² + 3t³) C₁ + (3t² - 3t³) C₂ + t³ P₁ =
//
//        = P₀ + 3(C₁ - P₀)t + 3(C₂ - 2C₁ + P₀)t² + (P₁ - 3C₂ + 3C₁ - P₀)t³
//
// First derivative
//
//   P'(t) = (-3 + 6t - 3t²) P₀ + (3 - 12t + 9t²) C₁ + (6t - 9t²) C₂ + 3t² P₁ =
//
//         = 3(C₁ - P₀) + 6(C₂ - 2C₁ + P₀) t + 3(P₁ - 3C₂ + 3C₁ - P₀) t²
//
class BezierCurve
{
 private:
  Vector P0_;           // Start point.
  Vector C1_;           // First control point.
  Vector C2_;           // Second control point.
  Vector P1_;           // End point.

 public:
  // Construct an incomplete BezierCurve with just the begin and end points set.
  BezierCurve(Vector P0, Vector P1) : P0_(P0), P1_(P1) { }

  // Construct a fully defined BezierCurve from begin and end point plus two control points.
  BezierCurve(Vector P0, Vector C1, Vector C2, Vector P1) : P0_(P0), C1_(C1), C2_(C2), P1_(P1) { }

  // Construct a fully defined BezierCurve from "matrix" columns.
  BezierCurve(BezierCurveMatrix const& m) :
    P0_(m.coefficient[0]),
    C1_(m.coefficient[1] / 3.0 + P0_),
    C2_(m.coefficient[2] / 3.0 + 2.0 * C1_ - P0_),
    P1_(m.coefficient[3] + 3.0 * (C2_ - C1_) + P0_) { }

  // Accessors.
  Vector P0() const { return P0_; }
  Vector C1() const { return C1_; }
  Vector C2() const { return C2_; }
  Vector P1() const { return P1_; }

  // As matrix (coefficients of a polynomial in t).
  BezierCurveMatrix M() const { return {P0_, 3.0 * (C1_ - P0_), 3.0 * (C2_ - 2.0 * C1_ + P0_), P1_ - 3.0 * (C2_ - C1_) - P0_}; }

  // Return the value at t.
  Vector P(double t) const
  {
    double t2 = t * t;
    double t3 = t2 * t;
    double one_minus_t = 1.0 - t;
    double one_minus_t_2 = one_minus_t * one_minus_t;
    double one_minus_t_3 = one_minus_t_2 * one_minus_t;

    return one_minus_t_3 * P0_ + 3.0 * one_minus_t_2 * t * C1_ + 3.0 * one_minus_t * t2 * C2_ + t3 * P1_;
  }

  // Velocity at t.
  Vector velocity(double t) const
  {
    double t2 = t * t;
    return 3.0 * ((-1 + 2 * t - t2) * P0_ + (1 - 4.0 * t + 3.0 * t2) * C1_ + (2.0 * t - 3.0 * t2) * C2_ + t2 * P1_);
  }

  // Special case: direction of the velocity at t=0.
  Direction D0() const
  {
    Vector velocity{3.0 * (C1_ - P0_)};
    return velocity.direction();
  }

  // Special case: direction of the velocity at t=1.
  Direction D1() const
  {
    Vector velocity{3.0 * (P1_ - C2_)};
    return velocity.direction();
  }

  // Calculate the axis aligned rectangle within which this BezierCurve exists (on the range 0 <= t <= 1).
  Rectangle extents() const;

  // Initialize a quadractic BezierCurve from two more data points.
  bool quadratic_from(Vector P_beta, Vector P_gamma);

  // Initialize a quadractic BezierCurve from one more data point and the direction vector in P₀.
  bool quadratic_from(Direction D0, Vector P_gamma);

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

} // namespace cairowindow
