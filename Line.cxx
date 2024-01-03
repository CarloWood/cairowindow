#include "sys.h"
#include "Line.h"
#include "debug.h"

namespace cairowindow {

Point Line::intersection_with(Line const& line2) const
{
  // Line1: P1 + λ D1
  // Line2: P2 + ξ D2

  Point const& P1 = point_;
  Direction const& D1 = direction_;
  Point const& P2 = line2.point();
  Direction const& D2 = line2.direction();

  // Let R2 be D2 rotated counter-clockwise by PI/2 (this is floating-point round off error free).
  Direction R2 = D2.normal();

  // Take dot product of D1 with R2:
  double D1_dot_R2 = D1.dot(R2);
  Dout(dc::notice, "D1_dot_R2 = " << D1_dot_R2);

  // If the lines are nearly parallel, use a different algorithm to calculate the intersection.
  if (std::abs(D1_dot_R2) < 0.01)
  {
    Direction d2{D2};

    // First make sure D1 and D2 point in (nearly) the same direction.
    if (D1.dot(D2) < 0.0)
      d2 = D2.inverse();

    // Let D = (D1 + D2)/2.            ^
    // Let R = (-Dy, Dx) / |D|         |
    //                                 |
    // Then,                           |
    //                                 |
    //   D1 = D + ε R1                 |R
    //   D2 = D - ε R2                 |
    //                                 |
    //                                 | D1
    //                                 |/
    //                          E ^    O                     D1+D2
    // --------O------------------|----+-----------------------O
    //       origin               '    O\
    //                                /  \
    //                              D2    D=(D1+D2)/2
    // Let E = D1 - D2
    double Ex = D1.x() - d2.x();
    double Ey = D1.y() - d2.y();
    //
//    double Dx = (D1.x() + D2.x()) / 2;
//    double Dy = (D1.y() + D2.y()) / 2;
//    Direction R(-Dy, Dx);

    // Calculate ε.
    double epsilon = sqrt(Ex * Ex + Ey * Ey) / 2;

    //Direction D12(D1.x() + d2.x(), D1.y() + d2.y());

    // Line1: P1 + λ D1 -->
    //
    //
  }

  //            intersection
  //                 \ /         P2
  //  --------+-------+-<--------+
  //       ^  |      /    D2  1  |
  //       |  |     /λ           |R2
  //     a |  |_  _/            1|
  //       | ^|   /|             |
  //       |b|| 1/               v
  //       | || /D1
  //       v v|/
  //          +P1
  //         /
  // intersection: P1 + λ D1 = P2 + ξ D2 -->
  //
  //   P2 - P1 = λ D1 - ξ D2
  //
  // Take dot product with R2 (D2 rotated PI/2 radians):
  //
  // λ = (P2 - P1)·R2 / (D1·R2)
  //
  // Note, in the above picture: a = -(P2 - P1)·R2, and b = -D1·R2

  // Take dot product of P2-P1 with R2:
  double P21_dot_R2 = (line2.point().x() - point_.x()) * R2.x() + (line2.point().y() - point_.y()) * R2.y();

  // Calculate lambda.
  double lambda = P21_dot_R2 / D1_dot_R2;

  // Return intersection point.
  return {point_.x() + lambda * direction_.x(), point_.y() + lambda * direction_.y()};
}

} // namespace cairowindow
