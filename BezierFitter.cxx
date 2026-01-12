#include "sys.h"
#include "BezierFitter.h"
#include "math/Vector.h"
#include <array>
#ifdef CWDEBUG
#include "utils/macros.h"
#endif

namespace cairowindow {

#ifdef CWDEBUG
std::string to_string(Orientation orientation)
{
  switch (orientation)
  {
    AI_CASE_RETURN(Orientation::horizontal);
    AI_CASE_RETURN(Orientation::vertical);
    AI_CASE_RETURN(Orientation::perpendicular);
  }
  AI_NEVER_REACHED
}
#endif

namespace {

constexpr double sqrt3 = 1.7320508075688773;
constexpr double c0 = sqrt3;
constexpr double d0 = 0.1666666666666666;       // 1/6
constexpr double d1 = 0.0833333333333333;       // 1/12

template<typename Vector>
void vals2coeffs7(std::array<Vector, 7>& out, Vector P0, Vector P1, Vector P2, Vector P3, Vector P4, Vector P5, Vector P6)
{
  Vector z0 = P6 + P0;
  Vector z1 = P6 - P0;
  Vector z2 = P5 + P1;
  Vector z3 = c0 * (P5 - P1);
  Vector z4 = P4 + P2;
  Vector z5 = P4 - P2;
  Vector z6 = 2.0 * P3;
  Vector w0 = z0 + z6;
  Vector w1 = z0 - z6;
  Vector w2 = z2 + z4;
  Vector w3 = z2 - z4;
  Vector w4 = z1 + z5;

  out[0] = d1 * (w0 + 2.0 * w2);
  out[1] = d0 * (w4 + z3);
  out[2] = d0 * (w1 + w3);
  out[3] = d0 * (z1 - 2.0 * z5);
  out[4] = d0 * (w0 - w2);
  out[5] = d0 * (w4 - z3);
  out[6] = d1 * (w1 - 2.0 * w3);
}

} // namespace

void BezierFitter::solve(std::function<Point(double)> const& func, IntersectRectangle<CS::plot> const& viewport, double tolerance,
    double const t0, double const t6, BezierCurve::Vector const P0, BezierCurve::Vector const P3, BezierCurve::Vector const P6)
{
  using Vector = BezierCurve::Vector;

  static constexpr double chebpt7_1 = 0.5 - 0.25 * sqrt3;
  static constexpr double chebpt7_5 = 0.5 + 0.25 * sqrt3;

  double td = t6 - t0;
  // Use .raw() here: we drop the knowledge about the fact that this is in CS::plot coordinates, and just work with math:: types.
  Vector P1{func(t0 + td * chebpt7_1).raw()};
  Vector P2{func(t0 + td * 0.25).raw()};
  Vector P4{func(t0 + td * 0.75).raw()};
  Vector P5{func(t0 + td * chebpt7_5).raw()};

  std::array<Vector, 7> coeffs;
  vals2coeffs7(coeffs, P0, P1, P2, P3, P4, P5, P6);

  double const xresid = std::abs(coeffs[4].x()) + std::abs(coeffs[5].x()) + std::abs(coeffs[6].x());    // To compare against the tolerance.
  double const yresid = std::abs(coeffs[4].y()) + std::abs(coeffs[5].y()) + std::abs(coeffs[6].y());
  double const xspread = std::abs(coeffs[1].x()) + std::abs(coeffs[2].x()) + std::abs(coeffs[3].x()) + xresid; // For finding a bounding box.
  double const yspread = std::abs(coeffs[1].y()) + std::abs(coeffs[2].y()) + std::abs(coeffs[3].y()) + yresid;

  // Don't abort because a meager six points are outside the graph.
  if (depth_ > 0)
  {
    // If the curve is entirely outside the viewport, don't bother rendering:
    if ((coeffs[0].x() + xspread < viewport.x1()) || (coeffs[0].x() - xspread > viewport.x2()) ||
        (coeffs[0].y() + yspread < viewport.y1()) || (coeffs[0].y() - yspread > viewport.y2()))
      return;

    // If the spread is very large in either dimension, we need
    // to be more careful about detecting the edges of the viewport.
    if (xspread > viewport.x2() - viewport.x1() || yspread > viewport.y2() - viewport.y1())
    {
      double const xmax = std::max({P0.x(), P1.x(), P2.x(), P3.x(), P4.x(), P5.x(), P6.x()});
      double const xmin = std::min({P0.x(), P1.x(), P2.x(), P3.x(), P4.x(), P5.x(), P6.x()});
      double const ymax = std::max({P0.y(), P1.y(), P2.y(), P3.y(), P4.y(), P5.y(), P6.y()});
      double const ymin = std::min({P0.y(), P1.y(), P2.y(), P3.y(), P4.y(), P5.y(), P6.y()});
      if (xmax < viewport.x1() || xmin > viewport.x2() || ymax < viewport.y1() || ymin > viewport.y2())
        return;
    }
  }

//  Dout(dc::notice, "xresid = " << xresid << "; yresid = " << yresid);

  // If we hit the desired tolerance, return a single bezier segment:
  if ((xresid < tolerance && yresid < tolerance) || depth_ == 10)
  {
//    Dout(dc::warning(depth_ == 10), "Reached max. depth!");

    // Alias degree 6 polynomial to degree 3:
    coeffs[0] += coeffs[6];
    coeffs[1] += coeffs[5];
    coeffs[2] += coeffs[4];

    // Convert from Chebyshev to Bernstein basis, and return:
    Vector Pt0 = coeffs[0] - 1.666666666666667 * coeffs[2];
    Vector Pt1 = 5.0 * coeffs[3] - 0.333333333333333 * coeffs[1];
    result_.emplace_back(P0, Pt0 + Pt1, Pt0 - Pt1, P6);

    return;
  }

  ++depth_;
  solve(func, viewport, tolerance, t0, 0.5 * (t0 + t6), P0, P2, P3);
  solve(func, viewport, tolerance, 0.5 * (t0 + t6), t6, P3, P4, P6);
  --depth_;
}

void BezierFitter::solve(std::function<Point(double)>&& func, Range const& domain, Rectangle const& viewport, double tolerance)
{
  using Vector = BezierCurve::Vector;

  // Clear result data, in case this object is being re-used.
  result_.clear();

  if (tolerance == -1.0)
    tolerance = 1e-5 * viewport.height();

  double t0 = domain.min();
  double t6 = domain.max();
  depth_ = 0;
  solve(func, viewport, tolerance, t0, t6, Vector{func(t0).raw()}, Vector{func(0.5 * (t0 + t6)).raw()}, Vector{func(t6).raw()});
}

void BezierFitter::solve(std::function<void(Point p, Vector v)> const& draw_line,
    std::function<Point(double)> const& P, std::function<Vector(double)> const& T,
    IntersectRectangle<CS::plot> const& viewport, double fraction, Orientation orientation,
    double t0, double t4, BezierCurve::Point P0, BezierCurve::Vector T0, BezierCurve::Point P2, BezierCurve::Point P4, BezierCurve::Vector T4)
{
  DoutEntering(dc::notice, "BezierFitter::solve(P, T, " << viewport << ", " << fraction << ", " << orientation << ", " <<
      t0 << ", " << t4 << ", " << P0 << ", " << T0 << ", " << P2 << ", " << P4 << ", " << T4 << ")");

  using Point = BezierCurve::Point;
  using Vector = BezierCurve::Vector;
  using Direction = BezierCurve::Direction;

  // 0 --1-- 2 --3-- 4
  // ↑       ↑       ↑
  // t0  ↑   |   ↑   t4
  double     t2         = 0.5 * (t0 + t4);  // 2 is in between 0 and 4.
  double t1             = 0.5 * (t0 + t2);  // 1 is in between 0 and 2.
  double         t3     = 0.5 * (t2 + t4);  // 3 is in between 2 and 4.

  // Use .raw() here: we drop the knowledge about the fact that this is in CS::plot coordinates, and just work with math:: types.
  Point P1{P(t1).raw()};        // The point at t1.
  Vector T2{T(t2).raw()};       // The tangent of P2 (P2 itself is provided as input).
  Point P3{P(t3).raw()};        // The point at t3.

  if (depth_ >= 1)
  {
    result_.emplace_back(P0, P4);
    Vector S = result_.back().cubic_from(T0, T4, P2);

    // draw_line needs to be called with cairowindow::Point / cairowindow::Vector again.
    cairowindow::Point const cs_P0{P0};
    cairowindow::Vector const cs_V0{result_.back().V0()};
    draw_line(cs_P0, cs_V0);

    Vector D1 =
      orientation == Orientation::horizontal ? Direction::right
                                             : orientation == Orientation::vertical ? Direction::up
                                                                                    : T(t1).raw().rotate_90_degrees();
    Vector D3 =
      orientation == Orientation::horizontal ? Direction::right
                                             : orientation == Orientation::vertical ? Direction::up
                                                                                    : T(t3).raw().rotate_90_degrees();
    double max_distance_squared =
      utils::square(fraction * (orientation == Orientation::vertical ? viewport.y2() - viewport.y1() : viewport.x2() - viewport.x1()));

    bool reject = S.x() <= 0.0 || S.y() <= 0.0 ||
      result_.back().distance_squared(P1, D1) > max_distance_squared ||
      result_.back().distance_squared(P3, D3) > max_distance_squared;

    if (!reject)
    {
#if 0
      auto& bc = result_.front();
      Dout(dc::notice, "Drawing: " << bc << "; P0 = " << bc.P0() << ", V0 = " << bc.V0() << ", T0 = " << bc.V0().direction());
      Dout(dc::notice, "Direction: " << (bc.B(0.001) - bc.P0()).direction());
#endif
      return;
    }

    result_.pop_back();
  }

  ++depth_;
  solve(draw_line, P, T, viewport, fraction, orientation, t0, t2, P0, T0, P1, P2, T2);
  solve(draw_line, P, T, viewport, fraction, orientation, t2, t4, P2, T2, P3, P4, T4);
  --depth_;
}

void BezierFitter::solve(std::function<void(Point p, Vector v)> const& draw_line,
    std::function<Point(double)>&& P, std::function<Vector(double)>&& T,
    Range const& domain, Rectangle const& viewport, double fraction, Orientation orientation)
{
  // Clear result data, in case this object is being re-used.
  result_.clear();

  double t0 = domain.min();
  double t1 = domain.max();
  depth_ = 0;
  // Use .raw() here: we drop the knowledge about the fact that this is in CS::plot coordinates, and just work with math:: types.
  solve(draw_line, P, T, viewport, fraction, orientation, t0, t1, P(t0).raw(), T(t0).raw(), P(0.5 * (t0 + t1)).raw(), P(t1).raw(), T(t1).raw());
}

} // namespace cairowindow
