#pragma once

#include "Vector.h"
#include "Rectangle.h"
#include "utils/has_print_on.h"
#include <array>
#include <memory>
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
//   P(t) = (1-t)³ P₀ + 3(1-t)²t C₀ + 3(1-t)t² C₁ + t³ P₁
//
//        = (1 - 3t + 3t² - t³) P₀ + (3t - 6t² + 3t³) C₀ + (3t² - 3t³) C₁ + t³ P₁ =
//
//        = P₀ + 3(C₀ - P₀)t + 3(C₁ - 2C₀ + P₀)t² + (P₁ - 3C₁ + 3C₀ - P₀)t³
//
// First derivative
//
//   P'(t) = (-3 + 6t - 3t²) P₀ + (3 - 12t + 9t²) C₀ + (6t - 9t²) C₁ + 3t² P₁ =
//
//         = 3(C₀ - P₀) + 6(C₁ - 2C₀ + P₀) t + 3(P₁ - 3C₁ + 3C₀ - P₀) t²
//
// Second derivative
//
//   P''(t) = 6(1 - t) P₀ + 6(-2 + 3t) C₀ + 6(1 - 3t) C₁ + 6t P₁ =
//
//          = 6(C₁ - 2C₀ + P₀) + 6(P₁ - 3C₁ + 3C₀ - P₀) t
//
class BezierCurve
{
 private:
  Vector P0_;           // Start point.
  Vector C0_;           // First control point.
  Vector C1_;           // Second control point.
  Vector P1_;           // End point.

 public:
  BezierCurve() = default;

  // Construct an incomplete BezierCurve with just the begin and end points set.
  BezierCurve(Point P0, Point P1) : P0_(P0), P1_(P1) { }
  BezierCurve(Vector P0, Vector P1) : P0_(P0), P1_(P1) { }

  // Construct a fully defined BezierCurve from begin and end point plus two control points.
  BezierCurve(Vector P0, Vector C0, Vector C1, Vector P1) : P0_(P0), C0_(C0), C1_(C1), P1_(P1) { }
  BezierCurve(Point P0, Point C0, Point C1, Point P1) : P0_(P0), C0_(C0), C1_(C1), P1_(P1) { }

  // Construct a fully defined BezierCurve from "matrix" columns.
  BezierCurve(BezierCurveMatrix const& m) :
    P0_(m.coefficient[0]),
    C0_(m.coefficient[1] / 3.0 + P0_),
    C1_(m.coefficient[2] / 3.0 + 2.0 * C0_ - P0_),
    P1_(m.coefficient[3] + 3.0 * (C1_ - C0_) + P0_) { }

  // Accessors.
  Vector P0() const { return P0_; }
  Vector C0() const { return C0_; }
  Vector C1() const { return C1_; }
  Vector P1() const { return P1_; }

  // As matrix (coefficients of a polynomial in t).
  BezierCurveMatrix M() const { return {P0_, V0(), 0.5 * A0(), 1.0 / 6.0 * J0()}; }

  // Velocity at t=0.
  Vector V0() const
  {
    return 3.0 * (C0_ - P0_);
  }

  // Acceleration at t=0.
  Vector A0() const
  {
    return 6.0 * ((C1_ - P0_) - 2.0 * (C0_ - P0_));
  }

  // Jolt at t=0.
  Vector J0() const
  {
    return 6.0 * ((P1_ - P0_) - 3.0 * (C1_ - C0_));
  }

  // Point at t.
  Point P(double t) const
  {
    Vector B{P0_ + t * (V0() + t / 2.0 * (A0() + t / 3.0 * J0()))};
    return B.point();
  }

  // Velocity at t.
  Vector velocity(double t) const
  {
    return V0() + t * (A0() + t / 2.0 * J0());
  }

  // Acceleration at t.
  Vector acceleration(double t) const
  {
    // 6(C₁ - 2C₀ + P₀) + 6(P₁ - 3C₁ + 3C₀ - P₀) t
    return A0() + t * J0();
  }

  // Special case: direction of the velocity at t=0.
  Direction D0() const
  {
    return V0().direction();
  }

  // Special case: direction of the velocity at t=1.
  Direction D1() const
  {
    Vector V1{3.0 * (P1_ - C1_)};
    return V1.direction();
  }

  // Curvature at t.
  Vector curvature(double t) const
  {
    Vector v = velocity(t);
    Vector a = acceleration(t);
    return a.dot(v.rotate_90_degrees()) / utils::square(v.length_squared()) * v.rotate_90_degrees();
  }

  double chord_length(double tolerance) const;

  // Calculate the axis aligned rectangle within which this BezierCurve exists (on the range 0 <= t <= 1).
  Rectangle extents() const;

  // Initialize a quadractic BezierCurve from two more data points.
  bool quadratic_from(Vector P_beta, Vector P_gamma);
  bool quadratic_from(Point P_beta, Point P_gamma) { return quadratic_from(Vector{P_beta}, Vector{P_gamma}); }

  // Initialize a quadractic BezierCurve from one more data point and the direction vector in P₀.
  bool quadratic_from(Direction D0, Vector P_gamma);
  bool quadratic_from(Direction D0, Point P_gamma) { return quadratic_from(D0, Vector{P_gamma}); }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

namespace draw {
class BezierCurve;
} // namespace draw

namespace plot {
class Plot;

class BezierCurve : public cairowindow::BezierCurve
{
 public:
  using cairowindow::BezierCurve::BezierCurve;
  BezierCurve(cairowindow::BezierCurve const& bezier_curve) : cairowindow::BezierCurve(bezier_curve) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::BezierCurve> draw_object_;
};

} // namespace plot
} // namespace cairowindow
