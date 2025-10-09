#pragma once

#include "math/Hyperplane.h"
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

using math::Hyperplane;

#if 0
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
#endif

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
  using VectorType = math::Vector<n, FloatType>;
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
  FloatType distance_of(Corner corner, Hyperplane<n, FloatType> const& plane) const
  {
    return plane.distance(C_[corner.index()]);
  }

  std::vector<VectorType> intersection_points(Hyperplane<n, FloatType> const& plane);

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
  Side(HyperBlock<FloatType, n> const& block, Hyperplane<n, FloatType> const& plane) : s_(number_of_corners)
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
std::vector<typename HyperBlock<FloatType, n>::VectorType> HyperBlock<FloatType, n>::intersection_points(Hyperplane<n, FloatType> const& plane)
{
  DoutEntering(dc::notice, "HyperBlock<" << type_info_of<FloatType>().demangled_name() << ", " << n << ">::intersection_points(" << plane << ")");
  std::vector<VectorType> intersections;

  utils::Vector<math::Sign, CornerIndex> side(number_of_corners);
  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
    side[ci] = plane.side(C_[ci]);

  Dout(dc::notice|continued_cf, "side = ");
  for (CornerIndex ci = C_.ibegin(); ci != C_.iend(); ++ci)
    Dout(dc::continued, (side[ci] == math::positive ? "+" : side[ci] == math::in_plane ? "0" : "-"));
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

  return intersections;
}

} // namespace intersections
