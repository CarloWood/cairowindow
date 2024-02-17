#include "sys.h"
#include "BezierFitter.h"
#include "Vector.h"
#include <array>

namespace cairowindow {

#ifdef CWDEBUG
void BezierCurve::print_on(std::ostream& os) const
{
  os << "{P0:" << P0_ << ", C0:" << C0_ << ", C1:" << C1_ << ", P1:" << P1_ << '}';
}
#endif

namespace {

constexpr double sqrt3 = 1.7320508075688773;
constexpr double c0 = sqrt3;
constexpr double d0 = 0.1666666666666666;       // 1/6
constexpr double d1 = 0.0833333333333333;       // 1/12

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

void BezierFitter::solve(std::function<Point(double)> const& func, IntersectRectangle const& viewport, double tolerance,
    double const t0, double const t6, Vector const P0, Vector const P3, Vector const P6)
{
  static constexpr double chebpt7_1 = 0.5 - 0.25 * sqrt3;
  static constexpr double chebpt7_5 = 0.5 + 0.25 * sqrt3;

  double td = t6 - t0;
  Vector P1{func(t0 + td * chebpt7_1)};
  Vector P2{func(t0 + td * 0.25)};
  Vector P4{func(t0 + td * 0.75)};
  Vector P5{func(t0 + td * chebpt7_5)};

  std::array<Vector, 7> coeffs;
  vals2coeffs7(coeffs, P0, P1, P2, P3, P4, P5, P6);

  double const xresid = std::abs(coeffs[4].x()) + std::abs(coeffs[5].x()) + std::abs(coeffs[6].x());    // To compare against the tolerance.
  double const yresid = std::abs(coeffs[4].y()) + std::abs(coeffs[5].y()) + std::abs(coeffs[6].y());
  double const xspread = std::abs(coeffs[1].x()) + std::abs(coeffs[2].x()) + std::abs(coeffs[3].x()) + xresid; // For finding a bounding box.
  double const yspread = std::abs(coeffs[1].y()) + std::abs(coeffs[2].y()) + std::abs(coeffs[3].y()) + yresid;

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

//  Dout(dc::notice, "xresid = " << xresid << "; yresid = " << yresid);

  // If we hit the desired tolerance, return a single bezier segment:
  if (xresid < tolerance && yresid < tolerance)
  {
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

#ifdef CWDEBUG
  ++depth_;
  if (depth_ > 1000)
  {
    if (depth_ == 1001)
      Dout(dc::warning, "Reached max. depth!");
    return;
  }
#endif
  solve(func, viewport, tolerance, t0, 0.5 * (t0 + t6), P0, P2, P3);
  solve(func, viewport, tolerance, 0.5 * (t0 + t6), t6, P3, P4, P6);
#ifdef CWDEBUG
  --depth_;
#endif
}

void BezierFitter::solve(std::function<Point(double)>&& func, Range const& domain, Rectangle const& viewport, double tolerance)
{
  double t0 = domain.min();
  double t6 = domain.max();
#ifdef CWDEBUG
  depth_ = 0;
#endif
  solve(func, viewport, tolerance, t0, t6, Vector{func(t0)}, Vector{func(0.5 * (t0 + t6))}, Vector{func(t6)});
}

} // amespace cairowindow
