#include "sys.h"
#include "BezierCurve.h"
#include "math/CubicPolynomial.h"
#include "utils/square.h"
#include "utils/almost_equal.h"
#include "utils/macros.h"
#include <Eigen/Dense>
#include <boost/math/quadrature/gauss_kronrod.hpp>

#include "QuickGraph.h"
#include "cwds/Restart.h"

namespace cairowindow {

Rectangle BezierCurve::extents() const
{
  Vector const P0{this->P0()};
  Vector const P1{this->P1()};
  std::vector<double> x_candidates = { P0.x(), P1.x() };
  std::vector<double> y_candidates = { P0.y(), P1.y() };

  // If the curve extents outside of the rectangle spanned by begin and end point
  // then the first derivative of x(t) and/or y(t) must be zero at the edge of
  // the extent.
  //
  // Since the first derivative is given by V0 + A0 t + J/2 t² we can find the
  // zeroes by solving two quadratic equations, for x and y respectively.
  //
  // Let `a t² + b t + c = 0`.
  Vector a = 0.5 * J();
  Vector b = A0();
  Vector c = V0();

  double discriminant_x = utils::square(b.x()) - 4.0 * a.x() * c.x();
  double discriminant_y = utils::square(b.y()) - 4.0 * a.y() * c.y();

  if (discriminant_x >= 0.0)
  {
    double root = std::sqrt(discriminant_x);
    double t_minus = (-b.x() - root) / (2.0 * a.x());
    if (t_minus >= 0.0 && t_minus <= 1.0)
    {
      Point P_minus = P(t_minus);
      x_candidates.push_back(P_minus.x());
    }
    if (discriminant_x > 0.0)
    {
      double t_plus = (-b.x() + root) / (2.0 * a.x());
      if (t_plus >= 0.0 && t_plus <= 1.0)
      {
        Point P_plus = P(t_plus);
        x_candidates.push_back(P_plus.x());
      }
    }
  }
  if (discriminant_y >= 0.0)
  {
    double root = std::sqrt(discriminant_y);
    double t_minus = (-b.y() - root) / (2.0 * a.y());
    if (t_minus >= 0.0 && t_minus <= 1.0)
    {
      Point P_minus = P(t_minus);
      y_candidates.push_back(P_minus.y());
    }
    if (discriminant_y > 0.0)
    {
      double t_plus = (-b.y() + root) / (2.0 * a.y());
      if (t_plus >= 0.0 && t_plus <= 1.0)
      {
        Point P_plus = P(t_plus);
        y_candidates.push_back(P_plus.y());
      }
    }
  }

  double x_min = x_candidates[0];
  double x_max = x_min;
  for (int i = 1; i < x_candidates.size(); ++i)
  {
    if (x_candidates[i] < x_min)
      x_min = x_candidates[i];
    if (x_candidates[i] > x_max)
      x_max = x_candidates[i];
  }
  double y_min = y_candidates[0];
  double y_max = y_min;
  for (int i = 1; i < y_candidates.size(); ++i)
  {
    if (y_candidates[i] < y_min)
      y_min = y_candidates[i];
    if (y_candidates[i] > y_max)
      y_max = y_candidates[i];
  }
  return {x_min, y_min, x_max - x_min, y_max - y_min};
}

void BezierCurve::linear_from()
{
  // m_.coefficient[0] was initialized with P₀.
  Vector P0(m_.coefficient[0]);
  // m_.coefficient[1] was initialized with P₁.
  Vector P1(m_.coefficient[1]);

  // Linear.
  m_.coefficient[1] = P1 - P0;
  m_.coefficient[2] = Vector{0.0, 0.0};
}

bool BezierCurve::quadratic_from(Direction D0, Direction D1)
{
  Vector const P0{m_.coefficient[0]};
  Vector const P1{m_.coefficient[1]};

  Direction N1 = D1.rotate_90_degrees();
  double D0_dot_N1 = D0.dot(N1);

  if (std::abs(D0_dot_N1) < 10e-9)
  {
    // D0 and D1 are (almost) the same.
    if (std::abs((P1 - P0).dot(N1)) < 10e-9 && D0.dot(D1) > 0.0)
    {
      // Just draw a straight line.
      linear_from();
      return true;
    }
    return false;
  }

  // See Line::intersection_with.
  double P1P0_dot_N1 = (P1.x() - P0.x()) * N1.x() + (P1.y() - P0.y()) * N1.y();
  double lambda = P1P0_dot_N1 / D0_dot_N1;

  // This would lead to a negative velocity (aka, start off in the opposite direction of D0).
  if (lambda <= 0.0)
    return false;

  Vector C = P0 + Vector{D0, lambda};
  Vector V0{2.0 * (C - P0)};

  // The intersection point of lines through D0 and D1 are the control point C.
  // C = P0 + Vector{D0, lambda};

  // P(t) = (1 - t)² P₀ + 2(1 - t)t C + t² P₁
  //      = (1 - 2t + t²) P₀ + (2t - 2t²) C + t² P₁
  //      = P₀ + (2C - 2P₀) t + (P₀ - 2C + P₁) t²
  // V₀ = 2 (C - P₀)
  //quadratic_from({D0, 2 * lambda});
  quadratic_from(V0);

  // Make sure the exit velocity matches D1.
  Vector V1 = m_.coefficient[1] + 2.0 * m_.coefficient[2];
  if (V1.dot(D1) < 0.0)
    return false;

  return true;
}

bool BezierCurve::quadratic_from(Vector P_gamma, Direction Di, point_nt point, double& gamma)
{
  DoutEntering(dc::notice, "BezierCurve::quadratic_from(" << P_gamma << ", " << Di << ", " << point << ")");

  Vector const P0{m_.coefficient[0]};
  Vector const P1{m_.coefficient[1]};

  Vector Q_gamma(P_gamma - P0);
  Vector S(P1 - P_gamma);

  // As usual, P(t) = P₀ + V₀ t + ½ A₀ t²               - note that P₀ = P(0).
  // And thus  V(t) = V₀ + A₀ t                         - first derivative.
  //
  // Let Q(t) = P(t) - P₀ = V₀ t + ½ A₀ t²
  //
  // Let Q₁ = Q(1) = V₀ + 0.5 A₀ --> A₀ = 2 (Q₁ - V₀) -->
  //
  //     Q(t) = V₀ t + (Q₁ - V₀) t²
  //          = V₀ t (1 - t) + Q₁ t²
  //
  // Let Qᵧ = Q(γ) = V₀ γ (1 - γ) + Q₁ γ² -->
  //     V₀ = (Qᵧ - Q₁ γ²) / (γ (1 - γ)) = (1 / (γ (1 - γ))) Qᵧ - (γ / (1 - γ)) Q₁
  //
  // Let S = P1 - Pᵧ = (P1 - P₀) - (Pᵧ - P₀) = Q₁ - Qᵧ
  //
  //     V₀ = (1 / (γ (1 - γ))) Qᵧ - (γ / (1 - γ)) Q₁ + (γ / (1 - γ)) Qᵧ - (γ / (1 - γ)) Qᵧ
  //        = (1 / (γ (1 - γ)) - (γ / (1 - γ))) Qᵧ - (γ / (1 - γ)) S =
  //        = ((1 + γ) / γ) Qᵧ - (γ / (1 - γ)) S
  //
  //     Vᵧ = V(γ) = V₀ + A₀ γ = V₀ + 2 (Q₁ - V₀) γ
  //        = V₀ (1 - 2 γ) + 2 γ Q₁ - 2 γ Qᵧ + 2 γ Qᵧ
  //        = (1 - 2 γ) ((1 + γ) / γ) Qᵧ - (1 - 2 γ) (γ / (1 - γ)) S + 2 γ S + 2 γ Qᵧ
  //        = ((1 - 2 γ) ((1 + γ) / γ) + 2 γ) Qᵧ + (2 γ - (1 - 2 γ) (γ / (1 - γ))) S
  //        = ((1 - γ) / γ) Qᵧ + (γ / (1 - γ)) S
  //
  //     V₁ = V(1) = V₀ + A₀ = V₀ + 2 (Q₁ - V₀)
  //        = -V₀ + 2 Q₁
  //        = -((1 + γ) / γ) Qᵧ + (γ / (1 - γ)) S + 2 Q₁ - 2 Qᵧ + 2 Qᵧ
  //        = (2 - (1 + γ) / γ) Qᵧ + ((γ / (1 - γ)) + 2) S
  //        = -((1 - γ) / γ) Qᵧ + ((2 - γ) / (1 - γ)) S

  // We want to find a and b such that
  //
  //   Di = a Qᵧ + b S                                                          (1)
  //
  // Take dot product of (1) with Qᵧ:
  //
  //   Di · Qᵧ = a (Qᵧ · Qᵧ) + b (S · Qᵧ)                                       (2)
  //
  // and multiply (2) with (S · S):
  //
  //   (Di · Qᵧ) (S · S) = a (Qᵧ · Qᵧ) (S · S) + b (S · Qᵧ) (S · S)             (3)
  //
  // Take dot product of (1) with S:
  //
  //   Di · S = a (Qᵧ · S) + b (S · S)                                          (4)
  //
  // and multiply (4) with (S · Qᵧ):
  //
  //   (Di · S) (S · Qᵧ) = a (Qᵧ · S) (S · Qᵧ) + b (S · S) (S · Qᵧ)             (5)
  //
  // Subtract (5) from (3):
  //
  //   (Di · Qᵧ) (S · S) - (Di · S) (S · Qᵧ) =  a ((Qᵧ · Qᵧ) (S · S) - (Qᵧ · S) (S · Qᵧ))
  //
  //       (Di · Qᵧ) (S · S) - (Di · S) (S · Qᵧ)
  //   a = ---------------------------------------------------
  //         (Qᵧ · Qᵧ) (S · S) - (Qᵧ · S)²
  //
  // Multiply (2) with (Qᵧ · S):
  //
  //   (Di · Qᵧ) (Qᵧ · S) = a (Qᵧ · Qᵧ) (Qᵧ · S) + b (S · Qᵧ) (Qᵧ · S)          (6)
  //
  // Multiply (4) with (Qᵧ · Qᵧ):
  //
  //   (Di · S) (Qᵧ · Qᵧ) = a (Qᵧ · S) (Qᵧ · Qᵧ) + b (S · S) (Qᵧ · Qᵧ)          (7)
  //
  // Subtract (6) from (7):
  //
  //   (Di · S) (Qᵧ · Qᵧ) - (Di · Qᵧ) (Qᵧ · S) = b ((S · S) (Qᵧ · Qᵧ) - (S · Qᵧ) (Qᵧ · S))
  //
  //       (Di · S) (Qᵧ · Qᵧ) - (Di · Qᵧ) (Qᵧ · S)
  //   b = ---------------------------------------
  //        (Qᵧ · Qᵧ) (S · S) - (Qᵧ · S)²
  //
  //            (Di · Qᵧ) (S · S) - (Di · S) (S · Qᵧ)
  //   a / b = ---------------------------------------
  //           (Di · S) (Qᵧ · Qᵧ) - (Di · Qᵧ) (Qᵧ · S)

  double Di_dot_Q_gamma = Vector{Di}.dot(Q_gamma);
  double Di_dot_S = Vector{Di}.dot(S);
  double S_dot_Q_gamma = S.dot(Q_gamma);
  double a_enumerator = Di_dot_Q_gamma * S.norm_squared() - Di_dot_S * S_dot_Q_gamma;
  double b_enumerator = Di_dot_S * Q_gamma.norm_squared() - Di_dot_Q_gamma * S_dot_Q_gamma;
  double a_div_b = a_enumerator / b_enumerator;

  switch (point)
  {
    case point0:
    {
      // Di is the tangent at point P₀.

      // The velocity V₀ (at P₀) is given by (see above)
      //
      //   V₀ = ((1 + γ) / γ) Qᵧ - (γ / (1 - γ)) S
      //
      // which means since Di has the same direction as V₀
      //
      //   a / b = ((1 + γ) / γ) / -(γ / (1 - γ)) = (γ² - 1) / γ² -->
      //
      //   γ = sqrt(1 / (1 - a / b))
      //

      if (a_div_b >= 1.0 || a_enumerator < 0.0)
        return false;

      gamma = 1.0 / std::sqrt(1.0 - a_div_b);

      Vector V0 = ((1.0 + gamma) / gamma) * Q_gamma - (gamma / (1.0 - gamma)) * S;
      quadratic_from(V0);

      return true;
    }
    case point_gamma:
    {
      // Di is the tangent at point Pᵧ.

      //
      // The velocity Vᵧ (at Pᵧ) is given by (see above)
      //
      //   Vᵧ = (1 - γ) / γ Qᵧ + γ / (1 - γ) S
      //
      // which means since Di has the same direction as Vᵧ
      //
      //   a / b = ((1 - γ) / γ) / (γ / (1 - γ)) = ((1 - γ) / γ)² -->
      //
      //   γ = 1 / (sqrt(a / b) + 1)
      //

      if (a_div_b < 0.0 || a_enumerator < 0.0)
        return false;

      gamma = 1.0 / (1.0 + std::sqrt(a_div_b));

      Vector V0 = ((1.0 + gamma) / gamma) * Q_gamma - (gamma / (1.0 - gamma)) * S;
      quadratic_from(V0);

      return true;
    }
    case point1:
    {
      // Di is the tangent at point P₁.

      // The velocity V₁ (at P₁) is given by (see above)
      //
      //   V₁ = -((1 - γ) / γ) Qᵧ + ((2 - γ) / (1 - γ)) S
      //
      // which means since Di has the same direction as V₁
      //
      //   a / b = -((1 - γ) / γ) / ((2 - γ) / (1 - γ)) = -(1 - γ)² / (γ (2 - γ))
      //
      //   γ² - 2γ + 1/(1 - a/b) = 0
      //
      //   γ = 1 +/- sqrt(1 - 1/(1 - a/b))
      //
      //

      if (a_div_b > 0.0 || a_enumerator > 0.0)
        return false;

      gamma = 1.0 - std::sqrt(1.0 - 1.0 / (1.0 - a_div_b));

      Vector V0 = ((1.0 + gamma) / gamma) * Q_gamma - (gamma / (1.0 - gamma)) * S;
      quadratic_from(V0);

      return true;
    }
    default:
      // Not implemented.
      ASSERT(false);
  }
  AI_NEVER_REACHED
}

bool BezierCurve::quadratic_from(Vector P_beta, Vector P_gamma)
{
  Vector const P0{m_.coefficient[0]};
  Vector const P1{m_.coefficient[1]};

  // The method used here is by Michael A. Lachance, Michael A. and Arthur J. Schwartz,
  // as published in the article "Four point parabolic interpolation" 1991.
  // See https://deepblue.lib.umich.edu/handle/2027.42/29347

  Eigen::Matrix3d R{
      { P0.x(),         P0.y(), 1.0 },
      { P_beta.x(), P_beta.y(), 1.0 },
      { P1.x(),         P1.y(), 1.0 }
    };

  Eigen::RowVector3d P_gamma_1(P_gamma.x(), P_gamma.y(), 1.0);

  auto q = P_gamma_1 * R.inverse();

  int qs_equal_to_one = 0;
  for (int i = 0; i < 3; ++i)
    if (utils::almost_equal(q(i), 1.0, 10e-9))
      ++qs_equal_to_one;
  if (qs_equal_to_one >= 2 || q(0) * q(1) * q(2) >= 0.0)
    return false;       // There is no parabola going through all four points.

  // For P_gamma to be after P1 (t > 1) q(1) must be larger than 0
  // and q(2) must be larger than 1 and the discriminant must be subtracted.
  if (q(1) < 0.0 || q(2) < 1.0)
    return false;

  double discriminant = utils::square(-2.0 * q(1) * q(2)) - 4.0 * q(1) * (1.0 - q(1)) * q(2) * (1.0 - q(2));
  double beta = (2.0 * q(1) * q(2) - std::sqrt(discriminant)) / (2.0 * q(1) * (1.0 - q(1)));

  // Only draw parabola's that have P_beta before P0.
  if (beta > 0.0)
    return false;

  // V is defined as a "Vandermonde" Matrix,
  //
  //       ⎡0  0  1⎤
  //   V = ⎜β² β  1⎟
  //       ⎣1  1  1⎦
  //
  // But we need and calculate ⅓V⁻¹ directly:
  Eigen::Matrix3d V_inverse{
    {         beta - 1.0,  1.0,       -beta },
    {  1.0 - beta * beta, -1.0, beta * beta },
    { beta * beta - beta,  0.0,         0.0 }
  };
  V_inverse /= beta * (beta - 1.0);

  // Calculate V⁻¹R.
  Eigen::Matrix3d W = V_inverse * R;

  // The paper uses (t² t  1) V⁻¹R, where-as we use M (1  t  t²)ᵀ,
  // fortunately, the translation between the two matrixes is
  // straightforward.
  double m01 = W(1, 0);
  double m11 = W(1, 1);
  double m02 = W(0, 0);
  double m12 = W(0, 1);

  // Apply to the matrix.
  m_.coefficient[1] = Vector{m01, m11};
  m_.coefficient[2] = Vector{m02, m12};

  return true;
}

bool BezierCurve::quadratic_from(double v0qa, double v1qa)
{
  constexpr double to_loops = 1.0 / (2.0 * M_PI);
  double v0qa_in_loops = to_loops * v0qa;
  double v1qa_in_loops = to_loops * v1qa;
  double foo0 = v0qa_in_loops - std::floor(v0qa_in_loops);
  double foo1 = v1qa_in_loops - std::floor(v1qa_in_loops);

  if (std::abs(foo0 - foo1) < 0.5)
    return false;

  Vector const P0{m_.coefficient[0]};
  Vector const P1{m_.coefficient[1]};

  // Work with a rotation and translation invariant basis.
  Vector Q1 = P1 - P0;
  Vector N1 = Q1.rotate_90_degrees();

  // Note that P0 + 0.5 * V0 is the control point (which is at the intersection of the tangent lines at P0 and P1).
  // To get the length of |0.5 * V0| we can simply use the sine law:
  double v0_div_q1 = 2.0 * std::sin(v1qa) / std::sin(v1qa - v0qa);
  // Now calculate V0 by also giving it the correct direction.
  Vector V0 = v0_div_q1 * std::cos(v0qa) * Q1 + v0_div_q1 * std::sin(v0qa) * N1;

  quadratic_from(V0);
  return true;
}

#if 0
void print(double dfdg, double g, std::function<double(double)>&& f)
{
  constexpr double epsilon = 1e-3;
  Dout(dc::notice, "dFdg = " << dfdg << " (should have been " << ((f(g + epsilon) - f(g - epsilon)) / (2 * epsilon))  << ").");
}

void print2(double d2fdg2, double g, std::function<double(double)>&& f)
{
  constexpr double epsilon = 1e-4;
  Dout(dc::notice, "d²F/dg² = " << d2fdg2 << " (should have been " <<
      ((f(g + epsilon) + f(g - epsilon) - 2 * f(g)) / (epsilon * epsilon))  << ").");
}
#endif

namespace {

class CubicBezierCurveImpl
{
 public:
  using Point = BezierCurve::Point;
  using Vector = BezierCurve::Vector;

 private:
  BezierCurve* parent_;
  Vector const T0;
  Vector const T1;
  Point const Pg;

  // Get P0 and P1 from the values passed to the constructor.
  Point const P0{parent_->M().coefficient[0].as_point()};
  Point const P1{parent_->M().coefficient[1].as_point()};

  // Calculate the determinant of the matrix formed by the tangent vectors T0 and T1.
  //
  //   | T_0x  T_1x |
  //   | T_0y  T_1y | = T_0x * T_1y - T_0y * T_1x
  //
  double const T0xT1 = T0.cross(T1);
  Vector const V1 = P1 - P0;
  Vector const Vg = Pg - P0;

  double const V1xT0 = V1.cross(T0);
  double const VgxT0 = Vg.cross(T0);
  double const V1xT1 = V1.cross(T1);
  double const VgxT1 = Vg.cross(T1);

  // None of these are a function of g.
  double const cghx = V1xT1 * T0.x() - V1xT0 * T1.x();
  double const cgx = VgxT1 * T0.x() - VgxT0 * T1.x();
  double const chx = V1xT0 * T1.x();
  double const cex = VgxT0 * T1.x();
  double const cghy = V1xT1 * T0.y() - V1xT0 * T1.y();
  double const cgy = VgxT1 * T0.y() - VgxT0 * T1.y();
  double const chy = V1xT0 * T1.y();
  double const cey = VgxT0 * T1.y();
  double const cx = -2 * T0xT1 * V1.x();
  double const cy = -2 * T0xT1 * V1.y();

  // The following variables are a function of g.
  double h;                     // g²(2g - 3)
  double g1mg;                  // g(1 - g)
  double g1mg2;                 // g²(1 - g)² = g1mg²
  double gh;                    // g³(2g - 3) = g·h
  double fx;                    // cghx·g·h + cgx·g + chx·h + cex
  double fy;                    // cghy·g·h + cgy·g + chy·h + cex

  double Jp_x;                  // x-coordinate of J⁺, where J⁺ = (T0 × T1)/6 J.
  double Jp_y;                  // y-coordinate of J⁺.

  // First derivatives, additionally needed to calculate dn2Jpdg.
  double dg1mg2dg;              // ∂g1mg2/∂g
  double dhdg;                  // ∂h/∂g
  double dghdg;                 // ∂(gh)/∂g
  double dfxdg;                 // ∂fx/∂g
  double dfydg;                 // ∂fy/∂g
  double g1mgm4;                // (g(1 - g))⁻⁴ = g1mg⁻⁴ = g1mg2⁻²
  double dJp_xdg;               // ∂Jp_x/∂g
  double dJp_ydg;               // ∂Jp_y/∂g
  double dn2Jpdg;               // ∂∥J⁺∥²/∂g

  // Additionally needed to calculate ddn2Jpdgdg (mostly second derivatives).
  double dg1mgm4dg;             // ∂g1mgm4/∂g
  double ddg1mg2dgdg;           // ∂dg1mg2dg/∂g
  double ddghdgdg;              // ∂dghdg/∂g
  double ddhdgdg;               // ∂dhdg/∂g
  double ddfxdgdg;              // ∂dfxdg/∂g
  double ddfydgdg;              // ∂dfydg/∂g
  double ddJp_xdgdg;            // ∂dJp_xdg/∂g
  double ddJp_ydgdg;            // ∂dJp_ydg/∂g
  double ddn2Jpdgdg;            // ∂dn2Jpdg/∂g

 private:
#ifdef CWDEBUG
  friend struct DebugAccessToCubicBezierCurveImpl;
#endif
  friend Vector BezierCurve::cubic_from(Vector T0, Vector T1, Point Pg);
  [[gnu::always_inline]] inline void calculate_derivatives(double g);
  [[gnu::always_inline]] inline Vector initialize(double g);
  [[gnu::always_inline]] double derivative() const { return dn2Jpdg; }
  [[gnu::always_inline]] double second_derivative() const { return ddn2Jpdgdg; }

 public:
  CubicBezierCurveImpl(BezierCurve* parent, Vector T0_, Vector T1_, Point Pg_) : parent_(parent), T0(T0_), T1(T1_), Pg(Pg_) { }
};

#ifdef CWDEBUG
struct DebugAccessToCubicBezierCurveImpl : public CubicBezierCurve
{
  CubicBezierCurveImpl impl_;

  DebugAccessToCubicBezierCurveImpl(Point P0, Point P1, Vector T0, Vector T1, Point Pg) : CubicBezierCurve(P0, P1), impl_(this, T0, T1, Pg) { }

  void initialize_from_g(double g) final { impl_.initialize(g); }
  double derivative(double g) final { impl_.calculate_derivatives(g); return impl_.derivative(); }
  double second_derivative(double g) final { impl_.calculate_derivatives(g); return impl_.second_derivative(); }
};
#endif

} // namespace

void CubicBezierCurveImpl::calculate_derivatives(double g)
{
  using utils::square;

  h = g * g * (2 * g - 3);
  g1mg = g * (1 - g);
  g1mg2 = square(g1mg);
  gh = g * h;
  fx = gh * cghx + g * cgx + h * chx + cex;
  fy = gh * cghy + g * cgy + h * chy + cey;

  // The coordinates of J⁺. Since J⁺ = (T0xT1 / 6) J, minimizing J⁺ means minimizing J.
  Jp_x = fx / g1mg2 + cx;
  Jp_y = fy / g1mg2 + cy;

  // Calculate ∂g1mg2/∂g = ∂(g²(1-g)²)/∂g = ∂(g²-2g³+g⁴)/∂g = 2g - 6g² + 4g³ = g(2 - 6g + 4g²).
  dg1mg2dg = g * (2 - 6 * g + 4 * square(g));
  // Calculate ∂h/∂g = ∂(g²(2g - 3))/∂g = ∂(2g³ - 3g²)/∂g = 6g² - 6g = -6g(1 - g).
  dhdg = -6 * g1mg;
  // Calculate ∂(gh)/∂g = ∂(2g⁴ - 3g³)/∂g = 8g³ - 9g² = g²(8g - 9).
  dghdg = square(g) * (8 * g - 9);
  // Calculate ∂fx/∂g = ∂(cghx·gh + cgx·g + chx·h + cex)/∂g.
  dfxdg = cghx * dghdg + cgx + chx * dhdg;
  // Calculate ∂fy/∂g = ∂(cghy·gh + cgy·g + chy·h + cey)/∂g.
  dfydg = cghy * dghdg + cgy + chy * dhdg;
  // Calculate g1mg2⁻² = g1mg⁻⁴.
  g1mgm4 = 1.0 / square(g1mg2);
  // Calculate ∂Jp_x/∂g = (∂fx/∂g · g1mg2 - fx · ∂g1mg2/∂g) / g1mg2².
  dJp_xdg = (dfxdg * g1mg2 - fx * dg1mg2dg) * g1mgm4;
  // Calculate ∂Jp_y/∂g = (∂fy/∂g · g1mg2 - fy · ∂g1mg2/∂g) / g1mg2².
  dJp_ydg = (dfydg * g1mg2 - fy * dg1mg2dg) * g1mgm4;
  // Calculate ∂∥J⁺∥²/∂g = ∂n2Jp/∂g = ∂(Jp_x² + Jp_y²)/∂g = ∂(Jp_x²)/∂g + ∂(Jp_y²)/∂g =
  //   2 Jp_x ∂Jp_x/∂g + 2 Jp_y ∂Jp_y/∂g
  dn2Jpdg = 2 * (Jp_x * dJp_xdg + Jp_y * dJp_ydg);

  // We need the derivative of dn2Jpdg (the second derivative of n2Jp). For that we need the derivative of each of the above factors.
  // Calculate ∂g1mgm4/∂g = ∂(g1mg2⁻²)/∂g = -2 g1mg2⁻³ dg1mg2dg.
  dg1mgm4dg = -2 * std::pow(g1mg2, -3.0) * dg1mg2dg;
  // Calculate ∂dg1mg2dg/∂g = ∂(2g - 6g² + 4g³)/∂g = 2 - 12g + 12g² = 2 - 12 g (1 - g).
  ddg1mg2dgdg = 2 - 12 * g1mg;
  // Calculate ∂dghdg/∂g = ∂(8g³ - 9g²)/∂g = 24g² - 18g = g (24g - 18).
  ddghdgdg = g * (24 * g - 18);
  // Calculate ∂dhdg/∂g = ∂(6g² - 6g)/∂g = 12g - 6.
  ddhdgdg = 12 * g - 6;
  // Calculate ∂dfxdg/∂g.
  ddfxdgdg = cghx * ddghdgdg + chx * ddhdgdg;
  // Calculate ∂dfydg/∂g.
  ddfydgdg = cghy * ddghdgdg + chy * ddhdgdg;
  // Calculate ∂dJp_xdg/∂g.
  ddJp_xdgdg = (ddfxdgdg * g1mg2 - fx * ddg1mg2dgdg) * g1mgm4 + (dfxdg * g1mg2 - fx * dg1mg2dg) * dg1mgm4dg;
  // Calculate ∂dJp_ydg/∂g.
  ddJp_ydgdg = (ddfydgdg * g1mg2 - fy * ddg1mg2dgdg) * g1mgm4 + (dfydg * g1mg2 - fy * dg1mg2dg) * dg1mgm4dg;
  // Calculate ∂dn2Jpdg/∂g.
  ddn2Jpdgdg = 2 * (square(dJp_xdg) + Jp_x * ddJp_xdgdg + square(dJp_ydg) + Jp_y * ddJp_ydgdg);

#if 0
  // Calculate the square of the norm of J⁺.
  double n2Jp = square(Jp_x) + square(Jp_y);

  Dout(dc::notice, "∥J⁺∥² = " << n2Jp << "; ∥J∥ = " << (6.0 * std::sqrt(n2Jp) / T0xT1));
  Dout(dc::notice, "∂(∥J⁺∥²)/∂g = " << dn2Jpdg << "; ∂²(∥J⁺∥²)/∂g² = " << ddn2Jpdgdg);
#endif
}

CubicBezierCurveImpl::Vector CubicBezierCurveImpl::initialize(double g)
{
  using utils::square;

  // We can call this function with any value of g (no need to first call calculate_derivatives).
  double const h = g * g * (2 * g - 3);

  // Helper vector.
  Vector const Z = h * V1 + Vg;

  //   Sx = (Z × T_1) / (3 (1-g)^2 g T0 × T1)
  //   Sy = (Z × T_0) / (3 g^2 (1-g) T0 × T1)
  Vector const S{
    Z.cross(T1) / (3 * square(1 - g) *      g  * T0xT1),
    Z.cross(T0) / (3 * square(g)     * (1 - g) * T0xT1)
  };

  // If not both elements of S are positive, then the Bezier curve enters P0 and/or P1 from the wrong end and this is not a usable solution.
  // Let R0 = C0 - P0 and R1 = C1 - P0; note that C1 - C0 = R1 - R0.
  Vector const R0 = S.x() * T0;
  Vector const R1 = V1 - S.y() * T1;

  auto& m_ = parent_->M();
  m_.coefficient[1] = 3.0 * R0;
  m_.coefficient[2] = 3.0 * (R1 - 2.0 * R0);
  m_.coefficient[3] = V1 - 3.0 * (R1 - R0);     // This is J/6.

  ASSERT(!std::isnan(m_.coefficient[0].x()));
  ASSERT(!std::isnan(m_.coefficient[0].y()));
  ASSERT(!std::isnan(m_.coefficient[1].x()));
  ASSERT(!std::isnan(m_.coefficient[1].y()));
  ASSERT(!std::isnan(m_.coefficient[2].x()));
  ASSERT(!std::isnan(m_.coefficient[2].y()));
  ASSERT(!std::isnan(m_.coefficient[3].x()));
  ASSERT(!std::isnan(m_.coefficient[3].y()));

  return S;
}

// Construct a Bezier curve defined by starting point P0 and end point P1 (initialized by the constructor),
// and the tangents T0 and T1;
//
//   B(t) = (1 - t)^3 * P0 + 3 * (1 - t)^2 * t * C0 + 3 * (1 - t) * t^2 * C1 + t^3 * P1
//
// where
//
//   C0 = P0 + Sx * T0;
//   C1 = P1 - Sy * T1;
//
// where Sx and Sy are two positive real values, such that B(g) = Pg.
//
// It is possible that either of the elements of the returned S is negative,
// in which case the curve is tangent to T0 and/or T1 respectively but
// approaching P0 and/or P1 from the wrong direction. The caller has to
// check this themselves.
//
BezierCurve::Vector BezierCurve::cubic_from(Vector T0, Vector T1, Point Pg)
{
  DoutEntering(dc::notice, "BezierCurve::cubic_from(" << T0 << ", " << T1 << ", " << Pg << ")");

  CubicBezierCurveImpl helper(this, T0, T1, Pg);

  // Iteratively change g until the change in g becomes less than one thousands of g.

  double prev_g;
  double g = 0.5;
  do
  {
    prev_g = g;

    helper.calculate_derivatives(g);
    double second_derivative = helper.second_derivative();

    // If the second derivative is negative, we'd go in the wrong direction (away from the root).
    bool use_bracket_algorithm = second_derivative <= 0.0;
    if (!use_bracket_algorithm)
    {
      // Newton-Raphson.
      g -= helper.derivative() / second_derivative;

      // If we went out of bounds, drop back to using a bracket algorithm.
      use_bracket_algorithm = g <= 0.0 || g >= 1.0;
    }
    if (use_bracket_algorithm)
    {
      double closest_boundary = prev_g < 0.5 ? 0.0 : 1.0;
      g = 0.5 * (closest_boundary + prev_g);
    }
//    Dout(dc::notice, (use_bracket_algorithm ? "Bracket" : "Newton-Raphson") <<
//        ": new g = " << g << "; Δg = " << (g - prev_g) << "; Δg/g = " << ((g - prev_g) / g));
  }
  while (std::abs(g - prev_g) > std::abs(0.001 * g));

  // Initialize the Bezier curve with the found value of g.
  return helper.initialize(g);
}

#ifdef CWDEBUG
//static
std::unique_ptr<CubicBezierCurve> CubicBezierCurve::create(Point P0, Point P1, Vector T0, Vector T1, Point Pg)
{
  return std::make_unique<DebugAccessToCubicBezierCurveImpl>(P0, P1, T0, T1, Pg);
}
#endif

double BezierCurve::quadratic_arc_length() const
{
  // Only call this for quadratic Bezier curves.
  ASSERT(m_.coefficient[3].x() == 0.0 && m_.coefficient[3].y() == 0.0);
  // Calculate the length with an algebraic formula for the quadratic Bezier.
  double v02 = V0().norm_squared();
  double v0 = std::sqrt(v02);
  double a02 = A0().norm_squared();
  double a0 = std::sqrt(a02);
  double z = V0().dot(A0());
  double s = std::sqrt(v02 + 2.0 * z + a02);
  double a03 = a02 * a0;
  return ((z * a0 + a03) * s - z * a0 * v0 + (v02 * a02 - z * z) * std::log((z + a02 + a0 * s) / (z + v0 * a0))) / (2.0 * a03);
}

double BezierCurve::quadratic_bending_energy() const
{
  // Only call this for quadratic Bezier curves.
  ASSERT(m_.coefficient[3].x() == 0.0 && m_.coefficient[3].y() == 0.0);

  // The bending energy is given by
  //        1
  //   E = ∫ 𝜅(t)² dt
  //      0
  // where 𝜅(t) is the curvature at t, which is given by
  //
  //           q'(t) x q"(t)
  //   𝜅(t) = ---------------
  //             |q'(t)|³
  //
  // where q(t) is the quadratic Bezier curve, q' is the first
  // derivative and q" the second derivative.
  //
  //   q(t)  = P₀ + V₀ t + ½ A₀ t²
  //   q'(t) =      V₀   +   A₀ t
  //   q"(t) =               A₀
  //
  // And thus q'(t) x q"(t) = (V₀ + A₀ t) x A₀ = V₀ x A₀.
  //
  //                   1       1
  //   E = (V₀ x A₀)² ∫ --------------- dt
  //                 0  (|V₀ + A₀ t|³)²
  //
  // Note that (|V₀ + A₀ t|³)² = (|V₀ + A₀ t|²)³ =
  //                           = ((v₀_x + a₀_x t)² + (v₀_y + a₀_y t)²)³
  //                           = ((a₀_x² + a₀_y²)t² + 2(v₀_x a₀_x + v₀_y a₀_y)t + (v₀_x² + v₀_y²))³
  //                           = (|A₀|² t² + 2 V₀·A₀ t + |V₀|²)³
  //
  double a02 = A0().norm_squared();
  double z = V0().dot(A0());
  double dz = 2.0 * z;
  double v02 = V0().norm_squared();
  double cs = utils::square(V0().cross(A0()));
  double a04 = utils::square(a02);
  double dz2 = utils::square(dz);
  double e = 4.0 * a02 * v02 - dz2;
  double e2 = utils::square(e);
  double se = std::sqrt(e);
  double i1 = (-(dz + 2.0 * a02) * (dz2 - 6.0 * a02 * dz - 2.0 * a02 * (5.0 * v02 + 3.0 * a02))) / utils::square(a02 + dz + v02) +
    24.0 * a04 * std::atan((dz + 2.0 * a02) / se) / se;
  double i0 = (-dz * (dz2 - 2.0 * a02 * 5.0 * v02)) / utils::square(v02) + 24.0 * a04 * std::atan(dz / se) / se;
  double result = cs * (i1 - i0) / (2.0 * e2);
  return result;
}

// A tolerance of less than 1e-15 probably won't work,
// due to round off errors of the double themselves.
double BezierCurve::arc_length(double tolerance) const
{
  double m00 = m_.coefficient[0].x();
  double m10 = m_.coefficient[0].y();
  double m01 = m_.coefficient[1].x();
  double m11 = m_.coefficient[1].y();
  double m02 = m_.coefficient[2].x();
  double m12 = m_.coefficient[2].y();
  double m03 = m_.coefficient[3].x();
  double m13 = m_.coefficient[3].y();
  double c0 = utils::square(m01) + utils::square(m11);
  double c1 = 4 * m01 * m02 + 4 * m11 * m12;
  double c2 = 4 * utils::square(m02) + 6 * m01 * m03 + 4 * utils::square(m12) + 6 * m11 * m13;
  double c3 = 12 * m02 * m03 + 12 * m12 * m13;
  double c4 = 9 * utils::square(m03) + 9 * utils::square(m13);
  auto f = [c0, c1, c2, c3, c4](double t){
    return std::sqrt(c0 + t * (c1 + t * (c2 + t * (c3 + t * c4))));
  };
  return boost::math::quadrature::gauss_kronrod<double, 31>::integrate(f, 0.0, 1.0, tolerance);
}

int BezierCurve::calculate_intersection_points(Point Q, Vector direction, std::array<double, 3>& intersection_point_t_values_out) const
{
  // Create a polynomial c0 + c1 x + c2 x^2 + c3 x^3.
  double c0 = m_.coefficient[0].dot(direction) - Vector{Q}.dot(direction);
  double c1 = m_.coefficient[1].dot(direction);
  double c2 = m_.coefficient[2].dot(direction);
  double c3 = m_.coefficient[3].dot(direction);
  math::CubicPolynomial<double> p(c0, c1, c2, c3);

  return p.get_roots(intersection_point_t_values_out);
}

double BezierCurve::distance_squared(Point Q, Vector direction) const
{
  // There are at most three intersection points.
  std::array<double, 3> intersection_point_t_values;
  int n = calculate_intersection_points(Q, direction, intersection_point_t_values);
  double min_distance_squared = std::numeric_limits<double>::infinity();

  if (AI_UNLIKELY(n == 0))
    return min_distance_squared;

  for (int i = 0; i < n; ++i)
  {
    double t = intersection_point_t_values[i];
    if (t <= 0.0 || t >= 1.0)
      continue;
    min_distance_squared = std::min(min_distance_squared, (Q - this->P(t)).norm_squared());
  }

  if (min_distance_squared == std::numeric_limits<double>::infinity())
  {
    Rectangle ext = extents();
    QuickGraph qg("No intersection points", "x", "y",
        {ext.offset_x(), ext.offset_x() + ext.width()}, {ext.offset_y(), ext.offset_y() + ext.width()});
    qg.add_function([&](double x) -> double { return ext.offset_y(); }, color::cornsilk);
    for (double t = 0.0; t <= 1.0; t += 0.001)
    {
      cairowindow::Point const p(this->P(t));
      qg.add_point(p, draw::PointStyle{}({.filled_shape = 15}));
    }
    cairowindow::Point const q(Q);
    qg.add_point(q, draw::PointStyle{}({.filled_shape = 5}));
    qg.wait_for_keypress();
  }

  // Returns infinity if no intersection point for a t in the range <0, 1> existed.
  return min_distance_squared;
}

double BezierCurve::stretching_energy(double tolerance) const
{
  return utils::square(arc_length(tolerance));
}

// Integral from 0 to 1 of the square of the curvature.
double BezierCurve::bending_energy(double tolerance) const
{
  auto f = [this](double t){
    return curvature(t).norm_squared();
  };
  return boost::math::quadrature::gauss_kronrod<double, 31>::integrate(f, 0.0, 1.0, tolerance);
}

#ifdef CWDEBUG
void BezierCurve::print_on(std::ostream& os) const
{
  os << "{m_x0:" << m_.coefficient[0] << ", m_x1:" << m_.coefficient[1] <<
    ", m_x2:" << m_.coefficient[2] << ", m_x3:" << m_.coefficient[3] << '}';
}
#endif

} // namespace cairowindow
