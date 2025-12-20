#pragma once

#include "LinePiece.h"
#include "math/Matrix.h"

namespace cairowindow::cs {

template<CS cs>
class Matrix;

template<CS cs>
struct MatrixTypes
{
  static constexpr int rows = 2;
  static constexpr int cols = 2;
  using scalar_type    = double;
  using vector_type    = Vector<cs>;
  using derived_type   = Matrix<cs>;

  template<int R, int C>
  requires (R == 2 && C == 2)
  using matrix_type    = Matrix<cs>;
};

template<CS cs>
class Matrix : public math::MatrixOps<MatrixTypes<cs>>
{
 private:
  math::Matrix<2, 2> raw_;

 public:
  // Construct an uninitialized Matrix.
  Matrix() = default;

  // Construct the Matrix `k I`.
  explicit MatrixData(double k) : raw_(k) { }

  // Explicit conversion from math::Matrix.
  explicit Matrix(math::Matrix<2, 2> const& raw) : raw_(raw) { }

  math::Matrix<2, 2>& raw() { return raw_; }
  math::Matrix<2, 2> const& raw() const { return raw_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << utils::to_string(cs) << ":" << raw_;
  }
#endif
};

} // namespace cairowindow::cs
