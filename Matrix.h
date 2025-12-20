#pragma once

#if 0 // FIXME: this can be removed since it is no longer used apparently.
#include "Vector.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace cairowindow {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

// A 2x2 matrix.
//
//   ⎡ m00_  m01_ ⎤
//   ⎣ m10_  m11_ ⎦
//
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

  // Construct a vector from two column vectors.
  Matrix(Vector const& col0, Vector const& col1) :
    m00_(col0.x()), m01_(col1.x()), m10_(col0.y()), m11_(col1.y()) { }

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

  Matrix transpose() const
  {
    return { m00_, m10_, m01_, m11_ };
  }

  friend Matrix operator*(Matrix const& lhs, Matrix const& rhs)
  {
    return {lhs.m00_ * rhs.m00_ + lhs.m01_ * rhs.m10_, lhs.m00_ * rhs.m01_ + lhs.m01_ * rhs.m11_,
            lhs.m10_ * rhs.m00_ + lhs.m11_ * rhs.m10_, lhs.m10_ * rhs.m01_ + lhs.m11_ * rhs.m11_};
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

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "[" << m00_ << ", " << m01_ << "; " << m10_ << ", " << m11_ << "]";
  }
#endif
};

} // namespace cairowindow
#endif
