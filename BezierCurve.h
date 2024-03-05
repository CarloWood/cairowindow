#pragma once

#include "Vector.h"
#include "Rectangle.h"
#include "utils/has_print_on.h"
#include "utils/square.h"
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

enum point_nt
{
  point_beta,
  point0,
  point1,
  point_gamma
};


// The cubic Bezier curve is the 2D parametric function:
//
//   P(t) = (1-t)³ P₀ + 3(1-t)²t C₀ + 3(1-t)t² C₁ + t³ P₁
//
//        = (1 - 3t + 3t² - t³) P₀ + (3t - 6t² + 3t³) C₀ + (3t² - 3t³) C₁ + t³ P₁ =
//
//        = P₀ + 3(C₀ - P₀)t + 3(C₁ - 2C₀ + P₀)t² + (P₁ - 3C₁ + 3C₀ - P₀)t³
//
//        = B  +        V0 t +          (A0/2) t² +                  J/6 t³
//
// The matrix m_ is defined as
//
//   m_ = [ B  V0  A0/2  J/6 ]
//
// First derivative
//
//   P'(t) = (-3 + 6t - 3t²) P₀ + (3 - 12t + 9t²) C₀ + (6t - 9t²) C₁ + 3t² P₁ =
//
//         = 3(C₀ - P₀) + 6(C₁ - 2C₀ + P₀) t + 3(P₁ - 3C₁ + 3C₀ - P₀) t²
//
//         = V0 + A0 t + J/2 t²
//
// Second derivative
//
//   P''(t) = 6(1 - t) P₀ + 6(-2 + 3t) C₀ + 6(1 - 3t) C₁ + 6t P₁ =
//
//          = 6(C₁ - 2C₀ + P₀) + 6(P₁ - 3C₁ + 3C₀ - P₀) t
//
//          = A0 + J t
//
class BezierCurve
{
 private:
  BezierCurveMatrix m_; // A 2x4 matrix where each colum represents a vector, respectively B, V0, A0/2 and J/6,
                        // where B is the position at t=0: P(0), V0 is the velocity at t=0: P'(0),
                        // A0 is the acceleration at t=0: P''(0) and J is the (constant) jolt: P'''(0).

 public:
  BezierCurve() = default;

  // Construct an incomplete BezierCurve with just the begin and end points set.
  // Put P₁ temporarily where V₀ will go.
  BezierCurve(Point P0, Point P1) : m_{{{Vector{P0}, Vector{P1}, {}, Vector{0.0, 0.0}}}} { }
  BezierCurve(Vector P0, Vector P1) : m_{{{P0, P1, {}, Vector{0.0, 0.0}}}} { }

  // Construct a fully defined BezierCurve from begin and end point plus two control points.
  BezierCurve(Point P0, Point C0, Point C1, Point P1) :
    m_{{{Vector{P0}, 3.0 * (C0 - P0), 3.0 * ((C1 - P0) - 2.0 * (C0 - P0)), P1 - P0 - 3.0 * (C1 - C0)}}} { }

  BezierCurve(Vector P0, Vector C0, Vector C1, Vector P1) :
    m_{{{P0, 3.0 * (C0 - P0), 3.0 * ((C1 - P0) - 2.0 * (C0 - P0)), P1 - P0 - 3.0 * (C1 - C0)}}} { }

  // Construct a fully defined BezierCurve from "matrix" columns.
  BezierCurve(BezierCurveMatrix const& m) : m_{m} { }

  // Accessors.
  Vector P0() const { return m_.coefficient[0]; }
  Vector C0() const { return m_.coefficient[1] / 3.0 + m_.coefficient[0]; }
  Vector C1() const { return m_.coefficient[2] / 3.0 + 2.0 / 3.0 * m_.coefficient[1] + m_.coefficient[0]; }
  Vector P1() const { return m_.coefficient[3] + m_.coefficient[2] + m_.coefficient[1] + m_.coefficient[0]; }

  // As matrix (coefficients of a polynomial in t).
  BezierCurveMatrix const& M() const { return m_; }

  // Velocity at t=0.
  Vector V0() const
  {
    return m_.coefficient[1];
  }

  // Acceleration at t=0.
  Vector A0() const
  {
    return 2.0 * m_.coefficient[2];
  }

  // Jolt at t=0.
  Vector J() const
  {
    return 6.0 * m_.coefficient[3];
  }

  // The vector to the starting point (B = Begin).
  Vector B(double t) const
  {
    return m_.coefficient[0] + t * (m_.coefficient[1] + t * (m_.coefficient[2] + t * m_.coefficient[3]));
  }

  // Point at t.
  Point P(double t) const
  {
    return B(t).point();
  }

  // Velocity at t.
  Vector velocity(double t) const
  {
    return V0() + t * (A0() + t * (0.5 * J()));
  }

  // Acceleration at t.
  Vector acceleration(double t) const
  {
    return A0() + t * J();
  }

  // Special case: direction of the velocity at t=0.
  Direction D0() const
  {
    return V0().direction();
  }

  // Special case: direction of the velocity at t=1.
  Direction D1() const
  {
    Vector V1{V0() + A0() + 0.5 * J()};
    return V1.direction();
  }

  // Curvature at t.
  Vector curvature(double t) const
  {
    Vector v = velocity(t);
    Vector a = acceleration(t);
    return a.dot(v.rotate_90_degrees()) / utils::square(v.length_squared()) * v.rotate_90_degrees();
  }

  double arc_length(double tolerance) const;
  double stretching_energy(double tolerance) const;
  double bending_energy(double tolerance) const;

  // Calculate the axis aligned rectangle within which this BezierCurve exists (on the range 0 <= t <= 1).
  Rectangle extents() const;

  // Initialize a "quadratic" BezierCurve (must use the (Point P0, Point P1) constructor) that represents a straight line between P₀ and P₁.
  void linear_from();

  // Initialize a quadratic BezierCurve from the velocity at P₀ (on top of going through the P₀ and P₁ passed to the constructor).
  void quadratic_from(Vector const& V0);

  // Initialize a quadratic BezierCurve from the direction of the tangents in P₀ and P₁.
  bool quadratic_from(Direction D0, Direction D1);

  // Initialize a quadratic BezierCurve from one more data point (Pᵧ) and the direction of the tangent in P_{point}.
  // Where point can be point0, point1 or point_gamma.
  bool quadratic_from(Vector P_gamma, Direction Di, point_nt point, double& gamma);
  bool quadratic_from(Point P_gamma, Direction Di, point_nt point, double& gamma) { return quadratic_from(Vector{P_gamma}, Di, point, gamma); }

  // Initialize a quadratic BezierCurve from two more data points (Pᵦ and Pᵧ).
  bool quadratic_from(Vector P_beta, Vector P_gamma);
  bool quadratic_from(Point P_beta, Point P_gamma) { return quadratic_from(Vector{P_beta}, Vector{P_gamma}); }

  // If this is a quadratic BezierCurve (constructed with quadratic_from) then this will return the exact arc length.
  // This uses an algebraic formula and therefore much faster than arc_length.
  double quadratic_arc_length() const;
  // Returns the square of the above.
  double quadratic_stretching_energy() const { return utils::square(quadratic_arc_length()); }
  // If this is a quadratic BezierCurve (constructed with quadratic_from) then this will return the exact bending energy.
  // This uses an algebraic formula and therefore much faster (and more accurate) than bending_energy.
  double quadratic_bending_energy() const;

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
