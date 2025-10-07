#pragma once

#include "utils/has_print_on.h"
#include "utils/ulong_to_base.h"
#include "utils/parity.h"
#include "utils/Vector.h"
#include "utils/macros.h"
#include "cwds/debug_ostream_operators.h"
#include <array>
#include <vector>
#include <concepts>
#include <stdexcept>
#include <iostream>
#include <iomanip>

namespace intersections {
// This class defines a print_on method.
using utils::has_print_on::operator<<;


// An n-dimensional vector (a column vector with n elements).
template<std::floating_point FloatType, int n>
struct Vector
{
  std::array<FloatType, n> v_;                          // The elements of the vector.

  // Construct a zeroed vector.
  Vector() : v_{} { }

  // Initialize this Vector from an initializer list.
  Vector(std::initializer_list<FloatType> v)
  {
    if (v.size() != n)
      throw std::invalid_argument("Initializer list must have exactly n elements.");
    std::copy(v.begin(), v.end(), v_.begin());
  }

  // Element access.
  FloatType& operator[](int i) { return v_[i]; }
  FloatType operator[](int i) const { return v_[i]; }

  // Add v2.
  Vector& operator+=(Vector const& v2)
  {
    for (int i = 0; i < n; ++i)
      v_[i] += v2[i];
    return *this;
  }

  // Add v1 and v2.
  friend Vector operator+(Vector const& v1, Vector const& v2)
  {
    Vector result(v1);
    result += v2;
    return result;
  }

  // Subtract v2.
  Vector& operator-=(Vector const& v2)
  {
    for (int i = 0; i < n; ++i)
      v_[i] -= v2[i];
    return *this;
  }

  // Subtract v2 from v1.
  friend Vector operator-(Vector const& v1, Vector const& v2)
  {
    Vector result(v1);
    result -= v2;
    return result;
  }

  // Return the dot product of v1 and v2.
  friend FloatType operator*(Vector const& v1, Vector const& v2)
  {
    FloatType result = 0;
    for (int i = 0; i < n; ++i)
      result += v1[i] * v2[i];
    return result;
  }

  // Elementwise multiply with scalar g.
  friend Vector operator*(FloatType g, Vector const& v)
  {
    Vector result(v);
    for (int i = 0; i < n; ++i)
      result[i] *= g;
    return result;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    char const* sep = "";
    os << '(';
    for (int i = 0; i < n; ++i)
    {
      os << sep << v_[i];
      sep = ", ";
    }
    os << ')';
  }
#endif
};

enum Sign
{
  negative = -1,      // The signed distance is negative: the point is on the -N side of the plane.
  in_plane = 0,       // The signed distance is very small: the point can be considered to be in the plane.
  positive = 1,       // The signed distance is positive: the point is on the +N side of the plane.
};

// Hyperplane
//
// An (n-1)-dimensional hyperplane orthogonal to a given normal vector N,
// located at the signed offset -b/‖N‖ from the origin.
//
// This class follows the convention originally set by Ludwig Otto Hesse, with respect to
// the signed distance. The Hesse form (https://en.wikipedia.org/wiki/Hesse_normal_form) of
// a hyperplane is where N is a unit normal vector (‖N‖ is one), indicated by giving N a hat,
// and reads:
//     ∧
//     N·X = ρ
//
// Here ρ is the signed offset of the plane from the origin: ρ is positive iff
// the plane lies in the direction that the normal points:
//
//              hyperplane
//                 |
//          N      |  N
//        O---->   |---->
//     origin      |
//                 P----------->A
//        -------->|  S(A) > 0
//           ρ > 0 |
//
// However, the signed distance S(A) from the plane to a point A is relative to the plane
// itself and positive if A is on the side of the plane that N is pointing towards
// (e.g. https://courses.csail.mit.edu/6.036/spring_2016/assignments/hw0_final.pdf (2d)).
//
// Let a point A be given by its projection P onto the plane and its signed distance s from
// the plane:
//               ∧
//     A = P + s N
//                                                 ∧
// Then the signed distance of A to the hyperplane N·X - ρ = 0 is:
//            ∧         ∧        ∧         ∧            ∧ ∧
//     S(A) = N·A - ρ = N·(P + s N) - ρ = (N·P - ρ) + s N·N = s
//
// because P is in the plane.
//
// Note that the signed distance of the origin itself is:
//            ∧
//     S(O) = N·O - ρ = -ρ,
//
// and has the opposite sign of the plane offset!
//
//              hyperplane
//                 |
//                 |  N
//     origin      |---->
//        O<-------Q-.         ∧
//            /    |  \_ Q = ρ·N
//   S(O) = -ρ < 0 |
//
// The constructor of this class is accepts the normal vector N (not necessarily
// a unit vector) as first argument and the negative of the dot product of N with
// some point P in the plane (b = -N·P).
//
// Note that in the literature it is common to write the equation of a hyperplane as:
//
//     N·X + b = 0,
//
// Hyperplane also stores N and b.
//
// Because P is in the plane, we have N·P + b = 0 --> b = -N·P.
// If we write the plane equation in Hesse form, by dividing by ‖N‖,
//
//     (N·X + b)/‖N‖ = 0 -> (N/‖N‖)·X = -b/‖N‖
//
// we see that b = -ρ ‖N‖, the plane offset multiplied by -‖N‖.
//
// Note:
// * The signed offset of the plane from the origin is ρ = -b/‖N‖.
//              ∧
// * The point ρN = (-b/‖N‖)(N/‖N‖) = (-b/‖N‖²)N is the projection of the origin onto the hyperplane.
//
// * The perpendicular distance between the origin and the plane is |ρ| = |b|/‖N‖.
//
// * The signed distance of the origin from the plane is -ρ = b/‖N‖.
//
template<std::floating_point FloatType, int n>
struct HyperPlane
{
  using VectorType = Vector<FloatType, n>;

  VectorType N_;                                        // The normal of the hyperplane.
  FloatType b_;                                         // The plane constant where N·X + b = 0.

  // Create a hyperplane that satisfies N·X + b = 0, where b = -N·P for some P on the plane.
  HyperPlane(VectorType const& N, FloatType b) : N_(N), b_(b) { }

  // Return the signed distance in Euclidean units (positive in +N direction).
  FloatType signed_distance(VectorType const& A) const
  {
    //        ∧
    // S(A) = N·A - ρ = (N·A + b) / ‖N‖.
    return (N_ * A + b_) / std::sqrt(N_ * N_);
  }

  // Return h such that A_projected = A - h N lies on the plane.
  FloatType height_along_N(VectorType const& A) const
  {
    // 0 = N·A_projected + b = N·(A - h N) + b = N·A + b - h N·N -->
    // h = (N·A + b) / N·N.
    return (N_ * A + b_) / (N_ * N_);
  }

  // Orthogonal projection of A onto the plane.
  VectorType project(VectorType const& A) const
  {
    return A - height_along_N(A) * N_;
  }

  // Return which side of the plane A is on.
  Sign side(VectorType const& A) const
  {
    FloatType h = height_along_N(A);
    FloatType ah = std::abs(h);
    constexpr FloatType abs_relative_error = 1e-6;
    if (AI_UNLIKELY(ah < abs_relative_error))
      return in_plane;
    return h > 0.0 ? positive : negative;
  }

  // Return intersection of the line through C1 and C2 with this HyperPlane.
  // Only call this function for points C1 and C2 where their difference is not perpendicular to N,
  // for example, where `distance` returns significantly different values for both points.
  VectorType intersection(VectorType const& C1, VectorType const& C2) const
  {
    // Let E be a line through C1 and C2: E: C1 + g(C2 - C1), where g parameterizes the points on E.
    // Fill that into the line equation to find the intersection:
    // N·(C1 + g(C2 - C1)) + b = 0 --> N·C1 + b + g N·(C2 - C1) = 0 --> g = -(N·C1 + b) / N·(C2 - C1)
    VectorType diff = C2 - C1;
    FloatType g = -(N_ * C1 + b_) / (N_ * diff);
    return C1 + g * diff;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{N:" << N_ << ", b:" << b_ << "}";
  }
#endif
};

class Corner;
using CornerIndex = utils::VectorIndex<Corner>;

// Don't use the most significant bit.
using value_type = uint32_t;

// A corner of a 31-cube.
// Also usable for n-cube' where n < 31: just leave all unused coordinates 0.
class Corner
{
 private:
  value_type c_;        // To be used as index into HyperBlock::C_ and Side::s_.

 public:
  // Construct the corner in the origin (0, 0, 0, ...).
  Corner() : c_{0} { }

  // Accessor.
  value_type value() const
  {
    return c_;
  }

  CornerIndex index() const
  {
    // Don't use the most significant bit.
    ASSERT((c_ & 0x80000000) == 0);
    return CornerIndex{c_};
  }

  void next()
  {
    // Advance to the next gray code.
    c_ ^= !utils::parity(c_) ? 1 : ((c_ & -c_) << 1);
  }

  Corner& operator++()
  {
    ++c_;
    return *this;
  }

  friend bool operator<(Corner lhs, Corner rhs)
  {
    return lhs.c_ < rhs.c_;
  }
};

template<std::floating_point FloatType, int n>
struct HyperBlock
{
  using VectorType = Vector<FloatType, n>;
  static constexpr int number_of_corners = 1 << n;

  utils::Vector<VectorType, CornerIndex> C_;          // The 2^n corners of the hyperblock.

  // Construct an axis-aligned hyperblock from two opposite corner vectors.
  HyperBlock(VectorType const& c1, VectorType const& c2) : C_(number_of_corners)
  {
    VectorType base = c2 - c1;
    for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
    {
      VectorType c;
      for (int d = 0; d < n; ++d)
      {
        int bit = 1 << d;
        if ((ci.get_value() & bit))
          c[d] += base[d];
      }
      C_[ci] = c1 + c;
    }
  }

  // Return the distance from the corner to the hyperplane.
  FloatType distance_of(Corner corner, HyperPlane<FloatType, n> const& plane) const
  {
    return plane.distance(C_[corner.index()]);
  }

  std::vector<VectorType> intersection_points(HyperPlane<FloatType, n> const& plane);

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    LIBCWD_USING_OSTREAM_PRELUDE
    os << "{C:" << C_ << "}";
  }
#endif
};

template<std::floating_point FloatType, int n>
class Side
{
 public:
  static constexpr int number_of_corners = HyperBlock<FloatType, n>::number_of_corners;

 private:
  utils::Vector<bool, CornerIndex> s_;        // Encodes which side a corner is on. The index is the corner.

 public:
  Side(HyperBlock<FloatType, n> const& block, HyperPlane<FloatType, n> const& plane) : s_(number_of_corners)
  {
    // Construct the origin corner.
    Corner current_corner;
    bool origin_distance_sign = block.distance_of(current_corner, plane) <= 0;
    // The origin corner is always on side 'false'.
    s_[current_corner.index()] = false;
    while (current_corner.index() == s_.iend())
    {
      ++current_corner;
      s_[current_corner.index()] = (block.distance_of(current_corner, plane) <= 0) != origin_distance_sign;
    }
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    for (CornerIndex ci = s_.ibegin(); ci != s_.iend(); ++ci)
      os << (s_[ci] ? '1' : '0');
  }
#endif
};

// Class that can print other objects.
class Printer
{
 private:
  int n_;       // The number of dimensions of the space that we live in.

 public:
  // Construct a Printer for an n-dimensional space.
  Printer(int n) : n_(n) { }

  void operator()(std::ostream& os, Corner const& corner) const
  {
    os << std::setw(n_) << std::setfill('0') << utils::ulong_to_base(corner.value(), "01");
  }
};

template<std::floating_point FloatType, int n>
std::vector<typename HyperBlock<FloatType, n>::VectorType> HyperBlock<FloatType, n>::intersection_points(HyperPlane<FloatType, n> const& plane)
{
  DoutEntering(dc::notice, "HyperBlock<" << type_info_of<FloatType>().demangled_name() << ", " << n << ">::intersection_points(" << plane << ")");
  std::vector<VectorType> intersections;

#if 1
  utils::Vector<Sign, CornerIndex> side(number_of_corners);
  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
    side[ci] = plane.side(C_[ci]);

  Dout(dc::notice|continued_cf, "side = ");
  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
    Dout(dc::continued, (side[ci] == positive ? "+" : side[ci] == in_plane ? "0" : "-"));
  Dout(dc::finish, "");

  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
  {
    for (int d = 0; d < n; ++d)
    {
      int bit = 1 << d;
      CornerIndex ci2(ci.get_value() | bit);
      if (side[ci] != side[ci2])
      {
        // Found two corners on opposite sides of the hyperplane.
        intersections.push_back(plane.intersection(C_[ci], C_[ci2]));
      }
    }
  }
#endif

  return intersections;
}

} // namespace intersections
