#include "sys.h"
#include "utils/square.h"
#include "BezierCurve.h"

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
      Vector P_minus = P(t_minus);
      x_candidates.push_back(P_minus.x());
    }
    if (discriminant_x > 0.0)
    {
      double t_plus = (-b.x() + root) / (2.0 * a.x());
      if (t_plus >= 0.0 && t_plus <= 1.0)
      {
        Vector P_plus = P(t_plus);
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
      Vector P_minus = P(t_minus);
      y_candidates.push_back(P_minus.y());
    }
    if (discriminant_y > 0.0)
    {
      double t_plus = (-b.y() + root) / (2.0 * a.y());
      if (t_plus >= 0.0 && t_plus <= 1.0)
      {
        Vector P_plus = P(t_plus);
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
  // Consider the quadratic Bezier curve:
  //
  //   P(t) = (1−t)² P₀ + 2t(1−t) C + t² P₁
  //
  // And thus
  //
  //   P(0) = P₀
  //   P(1) = P₁
  //
  // Furthermore, we have the two points Pᵦ and Pᵧ,
  // assumed to be part of the curve at t=β and t=γ respectively:
  //
  //   P(β) = (1−β)² P₀ + 2β(1−β) C + β² P₁ = Pᵦ
  //   P(γ) = (1−γ)² P₀ + 2γ(1−γ) C + γ² P₁ = Pᵧ
  //
  //   Pᵦ = (1 - 2β + β²) P₀ + (2β - 2β²) C + β² P₁ -->
  //   Pᵦ - P₀ = β (-(2 - β) P₀ + (2 - β) C - β C + β P₁) =
  //           = (2β - β²)(C - P₀) - β² ((C - P₀) - (P₁ - P₀)) =
  //           = 2(β - β²)(C - P₀) + β² (P₁ - P₀) -->
  //   C - P₀ = ((Pᵦ - P₀) - β² (P₁ - P₀)) / 2(β - β²)
  //
  // Lets apply a translation such that P₀ = (0, 0):
  // Let Pᵦ' = Pᵦ - P₀, P₁' = P₁ - P₀ and C' = C - P₀.
  // Then we have
  //
  //   C' = (Pᵦ' - β² P₁') / 2(β - β²)
  //
  // Likewise
  //
  //   C' = (Pᵧ' - γ² P₁') / 2(γ - γ²)
  //
  // Multiply both sides with 2(β - β²)(γ - γ²):
  //
  //   2(β - β²)(γ - γ²) C' = (γ - γ²)(Pᵦ' - β² P₁')
  //   2(β - β²)(γ - γ²) C' = (β - β²)(Pᵧ' - γ² P₁')
  //
  // Subtract
  //
  //   0 = (γ - γ²)Pᵦ' - (β - β²)Pᵧ' + (βγ² - γβ²)P₁'
  //
  // Let Nᵦ' be Pᵦ' rotated 90 degrees anti-clockwise.
  // Let Nᵧ' be Pᵧ' rotated 90 degrees anti-clockwise.
  // Let N₁' be P₁' rotated 90 degrees anti-clockwise.
  //
  // Take the dot product with Nᵦ' and Nᵧ' respectively.
  // Note that those must be independent or there won't
  // be a solution; aka, Nᵦ' and Nᵧ' form a basis.
  //
  // Since Pᵦ'·Nᵦ' = Pᵧ'·Nᵧ' = 0, we get the two equations:
  //
  //   0 = -(β - β²)(Pᵧ'·Nᵦ') + (βγ² - γβ²)(P₁'·Nᵦ')
  //   0 =  (γ - γ²)(Pᵦ'·Nᵧ') + (βγ² - γβ²)(P₁'·Nᵧ')
  //
  // Note that these dot products are just known constants;
  // lets define dᵧᵦ = Pᵧ'·Nᵦ', d₁ᵦ = P₁'·Nᵦ', dᵦᵧ = Pᵦ'·Nᵧ'
  // and d₁ᵧ = P₁'·Nᵧ'. So that we can abbreviate the above
  // equations as
  //
  //   0 = -(β - β²)dᵧᵦ + (βγ² - γβ²)d₁ᵦ
  //   0 =  (γ - γ²)dᵦᵧ + (βγ² - γβ²)d₁ᵧ
  //
  // Write the first one as polynomial in γ and the second
  // as a polynomial in β:
  //
  //   0 =  β d₁ᵦ γ² - β² d₁ᵦ γ - (β - β²) dᵧᵦ
  //   0 = -γ d₁ᵧ β² + γ² d₁ᵧ β + (γ - γ²) dᵦᵧ
  //
  // Lets divide the first one by β and the second one by -γ
  // (both are non-zero, or Pᵦ and/or Pᵧ would be equal to P₀).
  //
  //   0 = d₁ᵦ γ² - β d₁ᵦ γ + (β - 1) dᵧᵦ
  //   0 = d₁ᵧ β² - γ d₁ᵧ β + (γ - 1) dᵦᵧ
  //
  // Solving the latter:
  //
  //       γ d₁ᵧ +/- sqrt(γ² d₁ᵧ² - 4 d₁ᵧ (γ - 1) dᵦᵧ)
  //   β = -------------------------------------------
  //                      2 d₁ᵧ
  //
  // Note that β must be less than 0 by definition: Pᵦ
  // is supposed to become "before" P₀ (while Pᵧ comes
  // "after" P₁, aka γ > 1).
  //
  // Therefore the solution for β that is negative will be
  //
  //   β = (γ - sqrt(γ² - 4 (γ - 1) dᵦᵧ/d₁ᵧ))/2
  //
  // Let's abbreviate the square root for now with S, thus β = (γ - S)/2.
  // Injecting this into the former equation:
  //
  //   0 = d₁ᵦ γ² - ((γ - S)/2) d₁ᵦ γ + ((γ - S)/2 - 1) dᵧᵦ =
  //     = d₁ᵦ γ² - d₁ᵦ γ²/2 + (d₁ᵦ γ/2) S + (γ/2 - 1) dᵧᵦ - (dᵧᵦ/2) S =
  //     = d₁ᵦ γ²/2 + (γ/2 - 1) dᵧᵦ + (d₁ᵦ γ/2 - dᵧᵦ/2) S -->
  //
  // Multiply with two
  //
  //   0 = d₁ᵦ γ² + dᵧᵦ (γ - 2) + (d₁ᵦ γ - dᵧᵦ) S -->
  //
  //   S = (d₁ᵦ γ² + dᵧᵦ (γ - 2)) / (dᵧᵦ - d₁ᵦ γ)
  //
  // Squaring both sides
  //
  //   γ² - 4 (γ - 1) dᵦᵧ / d₁ᵧ = (d₁ᵦ γ² + dᵧᵦ (γ - 2))² / (dᵧᵦ - d₁ᵦ γ)²
  //

  return false;
}

bool BezierCurve::quadratic_from(Direction D0, Vector P_gamma)
{
  return false;
}

} // namespace cairowindow
