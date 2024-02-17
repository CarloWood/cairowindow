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

bool BezierCurve::quadratic_from(Vector P_beta, Vector P_gamma)
{
  Vector const P0{this->P0()};
  Vector const P1{this->P1()};

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
  m_.coefficient[3] -= m_.coefficient[2] + m_.coefficient[1];

  return true;
}

bool BezierCurve::quadratic_from(Direction D0, Vector P_gamma)
{
  DoutEntering(dc::notice, "BezierCurve::quadratic_from(" << D0 << ", " << P_gamma << ")");

  Vector const P0{this->P0()};
  Vector const P1{this->P1()};

  Vector P1P0(P0 - P1);
  Vector PgammaP0(P0 - P_gamma);

  // Test validity.
  {
    Vector PgammaP1(P1 - P_gamma);
    double P1R90dotPg = P1P0.rotate_90_degrees().dot(PgammaP0);
    double PgP1R90dotD0 = PgammaP1.rotate_90_degrees().dot(D0);

    if (P1R90dotPg > 0.0 || PgP1R90dotD0 > 0.0)
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

  double m00 = P0.x();
  double m10 = P0.y();
  double m01 = alpha * D0.x();
  double m11 = alpha * D0.y();
  double m02 = P1.x() - (m00 + m01);
  double m12 = P1.y() - (m10 + m11);

  // Apply to the matrix.
  m_.coefficient[1] = Vector{m01, m11};
  m_.coefficient[2] = Vector{m02, m12};
  m_.coefficient[3] -= m_.coefficient[2] + m_.coefficient[1];

  return true;
}

// A tolerance of less than 1e-15 probably won't work,
// due to round off errors of the double themselves.
double BezierCurve::chord_length(double tolerance) const
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

#ifdef CWDEBUG
void BezierCurve::print_on(std::ostream& os) const
{
  os << "{m_x0:" << m_.coefficient[0] << ", m_x1:" << m_.coefficient[1] <<
    ", m_x2:" << m_.coefficient[2] << ", m_x3:" << m_.coefficient[3] << '}';
}
#endif

} // namespace cairowindow
