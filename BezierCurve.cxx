#include "sys.h"
#include "utils/square.h"
#include "utils/almost_equal.h"
#include "BezierCurve.h"
#include <Eigen/Dense>

namespace cairowindow {

Rectangle BezierCurve::extents() const
{
  std::vector<double> x_candidates = { P0_.x(), P1_.x() };
  std::vector<double> y_candidates = { P0_.y(), P1_.y() };

  // If the curve extents outside of the rectangle spanned by begin and end point
  // then the first derivative of x(t) and/or y(t) must be zero at the edge of
  // the extent.
  //
  // Since the first derivative is given by 3(C₁ - P₀) + 6(C₂ - 2C₁ + P₀) t + 3(P₁ - 3C₂ + 3C₁ - P₀) t²
  // we can find the zeroes by solving two quadratic equations, for x and y respectively.
  //
  // Let `a t² + b t + c = 0`.

  Vector a = 3.0 * (P1_ - 3.0 * C2_ + 3.0 * C1_ - P0_);
  Vector b = 6.0 * (C2_ - 2.0 * C1_ + P0_);
  Vector c = 3.0 * (C1_ - P0_);

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

bool BezierCurve::quadratic_from(Vector P_beta, Vector P_gamma)
{
  // The method used here is by Michael A. Lachance, Michael A. and Arthur J. Schwartz,
  // as published in the article "Four point parabolic interpolation" 1991.
  // See https://deepblue.lib.umich.edu/handle/2027.42/29347

  Eigen::Matrix3d R(
      P0_.x(),       P0_.y(), 1.0,
      P_beta.x(), P_beta.y(), 1.0,
      P1_.x(),       P1_.y(), 1.0
    );

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
  Eigen::Matrix3d one_third_V_inverse(
              beta - 1.0,  1.0,       -beta,
       1.0 - beta * beta, -1.0, beta * beta,
      beta * beta - beta,  0.0,         0.0;
  one_third_V_inverse /= 3.0 * beta * (beta - 1.0);

  // Calculate ⅓V⁻¹R.
  Eigen::Matrix3d W = one_third_V_inverse * R;

  // The paper uses (t² t  1) V⁻¹R, where-as we use M (1  t  t²)ᵀ,
  // fortunately, the translation between the two matrixes is
  // straightforward. Note that we're still working with ⅓ of the
  // values because that is needed for the next formula.
  double one_third_m01 = W(1, 0);
  double one_third_m11 = W(1, 1);
  double one_third_m02 = W(0, 0);
  double one_third_m12 = W(0, 1);

  // Convert to control points.
  C1_ = Vector{one_third_m01, one_third_m11} + P0_;
  C2_ = Vector{one_third_m02, one_third_m12} + 2.0 * C1_ - P0_;

  return true;
}

bool BezierCurve::quadratic_from(Direction D0, Vector P_gamma)
{
  DoutEntering(dc::notice, "BezierCurve::quadratic_from(" << D0 << ", " << P_gamma << ")");

  Vector P1P0(P1_ - P0_);
  Vector PgammaP0(P_gamma - P0_);

  // Test validity.
  {
    Vector PgammaP1(P_gamma - P1_);
    double P1R90dotPg = P1P0.rotate_90_degrees().dot(PgammaP0);
    double PgP1R90dotD0 = PgammaP1.rotate_90_degrees().dot(D0);

    if (P1R90dotPg < 0.0 && PgP1R90dotD0 < 0.0)
      return false;
  }

  double PgR90dotD0 = PgammaP0.rotate_90_degrees().dot(D0);
  double P1R90dotD0 = P1P0.rotate_90_degrees().dot(D0);
  double gamma_squared = PgR90dotD0 / P1R90dotD0;
  double gamma = std::sqrt(gamma_squared);  // γ > 1 (comes after P₁).
  P1P0 *= gamma_squared;
  double PgdotP1 = PgammaP0.dot(P1P0);
  double alpha_squared =
    (PgammaP0.length_squared() + P1P0.length_squared() - 2 * PgdotP1) / (gamma_squared * utils::square(1.0 - gamma));
  double alpha = std::sqrt(alpha_squared);  // α > 0 (or we'd go in the wrong direction!)

  double m00 = P0_.x();
  double m10 = P0_.y();
  double m01 = alpha * D0.x();
  double m11 = alpha * D0.y();
  double m02 = P1_.x() - (m00 + m01);
  double m12 = P1_.y() - (m10 + m11);

  C1_ = Vector{m01, m11} / 3.0 + P0_;
  C2_ = Vector{m02, m12} / 3.0 + 2.0 * C1_ - P0_;

  return true;
}

} // namespace cairowindow
