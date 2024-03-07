#include "sys.h"
#include "utils/square.h"
#include "utils/almost_equal.h"
#include "BezierCurve.h"
#include <Eigen/Dense>
#include <boost/math/quadrature/gauss_kronrod.hpp>

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

  Direction N1 = D1.normal();
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
  double a_enumerator = Di_dot_Q_gamma * S.length_squared() - Di_dot_S * S_dot_Q_gamma;
  double b_enumerator = Di_dot_S * Q_gamma.length_squared() - Di_dot_Q_gamma * S_dot_Q_gamma;
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

double BezierCurve::quadratic_arc_length() const
{
  // Only call this for quadratic Bezier curves.
  ASSERT(m_.coefficient[3].x() == 0.0 && m_.coefficient[3].y() == 0.0);
  // Calculate the length with an algebraic formula for the quadratic Bezier.
  double v02 = V0().length_squared();
  double v0 = std::sqrt(v02);
  double a02 = A0().length_squared();
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
  double a02 = A0().length_squared();
  double z = V0().dot(A0());
  double dz = 2.0 * z;
  double v02 = V0().length_squared();
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

double BezierCurve::stretching_energy(double tolerance) const
{
  return utils::square(arc_length(tolerance));
}

// Integral from 0 to 1 of the square of the curvature.
double BezierCurve::bending_energy(double tolerance) const
{
  auto f = [this](double t){
    return curvature(t).length_squared();
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
