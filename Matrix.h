#pragma once

#include "Vector.h"

namespace cairowindow {

// A 2x2 matrix.
class Matrix
{
 private:
  double m00_;
  double m01_;
  double m10_;
  double m11_;

 public:
  Matrix(double m00, double m01, double m10, double m11) :
    m00_(m00), m01_(m01), m10_(m10), m11_(m11) { }

  Vector operator*(Vector const& v) const
  {
    return {m00_ * v.x() + m01_ * v.y(), m10_ * v.x() + m11_ * v.y()};
  }

  double determinant() const
  {
    return m00_ * m11_ - m10_ * m01_;
  }

  Matrix inverse() const
  {
    double const det = determinant();
    return { m11_ / det, -m10_ / det, -m01_ / det, m00_ / det };
  }

  friend Matrix operator*(double s, Matrix const& m)
  {
    Matrix result{m};
    result.m00_ *= s;
    result.m01_ *= s;
    result.m10_ *= s;
    result.m11_ *= s;
    return result;
  }
};

} // namespace cairowindow
