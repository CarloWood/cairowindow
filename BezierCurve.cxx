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
  Eigen::Matrix3d R;
  R << P0_.x(), P0_.y(), 1.0,
       P_beta.x(), P_beta.y(), 1.0,
       P1_.x(), P1_.y(), 1.0;

  Eigen::RowVector3d P_gamma_1;
  P_gamma_1 << P_gamma.x(), P_gamma.y(), 1.0;

  auto q = P_gamma_1 * R.inverse();

  int qs_equal_to_one = 0;
  for (int i = 0; i < 3; ++i)
  {
    if (utils::almost_equal(q(i), 1.0, 10e-9))
      ++qs_equal_to_one;
  }
  if (qs_equal_to_one >= 2 || q(0) * q(1) * q(2) >= 0.0)
  {
    // There is no parabola going through all four points.
    return false;
  }

  // For P_gamma to be after P1 (t > 1) q(1) must be larger than 0 and q(2) must be larger than 1 and the discriminant must be subtracted.
  if (q(1) < 0.0 || q(2) < 1.0)
    return false;

  double discriminant = utils::square(-2.0 * q(1) * q(2)) - 4.0 * q(1) * (1.0 - q(1)) * q(2) * (1.0 - q(2));
  double beta = (2.0 * q(1) * q(2) - std::sqrt(discriminant)) / (2.0 * q(1) * (1.0 - q(1)));

  // Only draw parabola's that have P1 before Q0 and P0.
  if (beta > 0.0)
    return false;

  Eigen::Matrix3d one_third_V_inverse;
  one_third_V_inverse << beta - 1.0, 1.0, -beta, 1.0 - beta * beta, -1.0, beta * beta, beta * beta - beta, 0.0, 0.0;
  one_third_V_inverse /= 3.0 * beta * (beta - 1.0);

  Eigen::Matrix3d W = one_third_V_inverse * R;

  double m01 = W(1, 0);
  double m11 = W(1, 1);
  double m02 = W(0, 0);
  double m12 = W(0, 1);

  C1_ = Vector{m01, m11} + P0_;
  C2_ = Vector{m02, m12} + 2.0 * C1_ - P0_;

  return true;
}

bool BezierCurve::quadratic_from(Direction D0, Vector P_gamma)
{
  return false;
}

} // namespace cairowindow
