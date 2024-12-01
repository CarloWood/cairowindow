#include "sys.h"
#include "cairowindow/QuickGraph.h"
#include "utils/square.h"
#include "utils/macros.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <complex>
#include <cassert>
#include <random>
#include <array>
#include <bitset>
#include "debug.h"

#if 0
void roots(std::complex<double> r0, std::complex<double> r1, std::complex<double> r2, cairowindow::QuickGraph& graph, std::array<cairowindow::draw::PointStyle, 35> const& point_styles)
{
  double c2 = std::real(-(r0 + r1 + r2));
  double c1 = std::real(r0 * r1 + r0 * r2 + r1 * r2);
  double c0 = std::real(-r0 * r1 * r2);
  double d = utils::square(c2) - 3.0 * c1;

//  double f2m30 = std::sqrt(std::abs(c0 * c0 / (c1 * c1 * c1)));

  if (d >= 0.0)
  {
    double sqrt_d = std::sqrt(d);
    double x_max = (-c2 - sqrt_d) / 3.0;
    double x_min = (-c2 + sqrt_d) / 3.0;

    double rr0 = std::real(r0);
    double rr1 = std::real(r1);
    double rr2 = std::real(r2);

    if (x_min < 0.0)
      ASSERT(rr2 < rr1 && rr1 < rr0);
    else if (0.0 < x_max)
      ASSERT(rr0 < rr1 && rr1 < rr2);
    else if (c2 < 0)
      ASSERT(rr1 < rr0 && rr0 < rr2);
    else
      ASSERT(rr2 < rr0 && rr0 < rr1);

//  if (std::isfinite(x) && -8.0 < x && x < 4.0)
//    graph.add_point({x, y}, point_styles[0]);
  }
}
#else

enum Order
{
  A,
  B,
  C,
  D
};

std::string to_string(Order order)
{
  switch (order)
  {
    AI_CASE_RETURN(A);
    AI_CASE_RETURN(B);
    AI_CASE_RETURN(C);
    AI_CASE_RETURN(D);
  }
  AI_NEVER_REACHED
}

std::ostream& operator<<(std::ostream& os, Order order)
{
  return os << to_string(order);
}

constexpr int A_bit = 1;
constexpr int B_bit = 2;
constexpr int C_bit = 4;
constexpr int D_bit = 8;
std::array<int, 64> table;

void roots(double r0, double r1, double r2)
{
  ASSERT(std::abs(r0) <= std::abs(r1) && std::abs(r1) <= std::abs(r2));

  // Coefficients of the monic cubic polynomial.
  double c2 = -(r0 + r1 + r2);
  double c1 = r0 * r1 + r0 * r2 + r1 * r2;
  double c0 = -r0 * r1 * r2;

  // Discriminant of the derivative
  double d = c2 * c2 - 3.0 * c1;
  ASSERT(d >= 0.0);

  double sqrt_d = std::sqrt(d);

  // Correct calculation of critical points
  double x_max = (-c2 - sqrt_d) / 3.0;
  double x_min = (-c2 + sqrt_d) / 3.0;

  // Determine which root to approximate.
  double ratio = 1.0;   // |root|/|r0|
  if (x_max > 0.0)
  {
    // Calculate the root that is less than x_max.
    int m = (r0 < x_max ? 1 : 0) | (r1 < x_max ? 2 : 0) | (r2 < x_max ? 4 : 0);
    ASSERT(m == 1);
  }
  else if (x_min < 0.0)
  {
    // Calculate the root that is larger than x_min.
    int m = (r0 > x_min ? 1 : 0) | (r1 > x_min ? 2 : 0) | (r2 > x_min ? 4 : 0);
    ASSERT(m == 1);
  }
  else
  {
    // Calculate the root that is larger than x_max but less than x_min.
    int m = (x_max < r0 && r0 < x_min ? 1 : 0) | (x_max < r1 && r1 < x_min ? 2 : 0) | (x_max < r2 && r2 < x_min ? 4 : 0);
    ASSERT(m == 1 || m == 2);
    if (m == 2)
      ratio = std::abs(r1) / std::abs(r0);
  }

  static double worst_case = 0.0;
  if (ratio > worst_case)
  {
    worst_case = ratio;
    Dout(dc::notice, worst_case << " : " << r0 << ", " << r1 << ", " << r2);
  }
}
#endif

// b = r₀⋅r₁ + r₀⋅r₂ + r₁⋅r₂
// c = - r₀⋅r₁⋅r₂
int main()
{
  Debug(NAMESPACE_DEBUG::init());

  std::streamsize const old_precision = std::cout.precision(15);
  std::ios::fmtflags const old_flags = std::cout.setf(std::ios::fixed);

  std::vector<double> values;

//  for (int factor4 = -3; factor4 <= 3; ++factor4)
    for (int factor3 = -3; factor3 <= 3; ++factor3)
      for (int factor2 = -3; factor2 <= 3; ++factor2)
        for (int factor1 = -3; factor1 <= 3; ++factor1)
          for (int factor0 = -3; factor0 <= 3; ++factor0)
          {
            double val = /*factor4 * 10000.0 +*/ factor3 * 1000.0 + factor2 * 100.0 + factor1 * 10.0 + factor0;
            values.push_back(val);
          }
  std::sort(values.begin(), values.end(), [](double a, double b){ return std::abs(a) < std::abs(b); });
  int s = values.size();

  //roots(3.0, -9.0, -3168.0);
  //return 0;

  cairowindow::QuickGraph graph("Small root", "log10(sqrt(|c0^2 / c1^3|))", "log10(|r0|/|r1|)", {-4.0, 4.0}, {-4.0, 0.0});
  graph.add_function([](double x){ return x; }, cairowindow::color::cornsilk);

  cairowindow::draw::PointStyle point_style_a({.color_index = 0, .filled_shape = 15});
  cairowindow::draw::PointStyle point_style_b({.color_index = 2, .filled_shape = 15});
  cairowindow::draw::PointStyle point_style_c({.color_index = 3, .filled_shape = 15});

  std::array<cairowindow::draw::PointStyle, 35> point_styles;
  for (int i = 0; i < point_styles.size(); ++i)
    point_styles[i] = cairowindow::draw::PointStyle({.color_index = i % 32, .filled_shape = 15});

  // Let P(x) = x³ - (r₀ + r₁ + r₂)⋅x² + (r₀⋅r₁ + r₀⋅r₂ + r₁⋅r₂)⋅x - r₀⋅r₁⋅r₂
  // where |r₀| ⩽ |r₁| ⩽ |r₂|.
  //
  // If the cubic has only one real roots, then the other two a complex conjugate
  // with equal absolute values. Hence, the following possibilities exist:
  //
  // a) P has three real roots (where |r₀| ⩽ |r₁| ⩽ |r₂|).
  // b) only r₀ is real and |r₀| ⩽ |r₁| = |r₂|
  // c) only r₂ is real and |r₀| = |r₁| ⩽ |r₂|
  // d) only r₁ is real and |r₀| = |r₁| = |r₂|, but in this case, without loss of generality, we swap r₁ and r₂ so that this becomes case c).
  //
  int seed = 10248;
  Dout(dc::notice, "seed = " << seed);
  std::mt19937 engine(seed);

  std::uniform_real_distribution<double> dist_less_than_ten(0.0001, 10.0);

  for (int case_abc = 0; case_abc < 1; ++case_abc)
  {
    switch (case_abc)
    {
      case 0:   // a
      {
        for (int i = 0; i < 100000000; ++i)
        {
//          if (i % 10000 == 0)
//            std::cout << i << '\n';
          double ar2 = dist_less_than_ten(engine);
          std::uniform_real_distribution<double> dist1(0.0001, ar2);
          double ar1 = dist1(engine);
          std::uniform_real_distribution<double> dist0(0.0001, ar1);
          double ar0 = dist0(engine);
          for (int signs = 0; signs <= 7; ++signs)
          {
            double r0 = (signs & 1) ? -ar0 : ar0;
            double r1 = (signs & 2) ? -ar1 : ar1;
            double r2 = (signs & 4) ? -ar2 : ar2;
            roots(r0, r1, r2 /*, graph, point_styles*/);
          }
        }
        break;
      }
#if 0
      case 1:   // b
        // r₀ is real, r₁ = α - 𝕚 β, r₂ = α + 𝕚 β.
        for (int i = 0; i < 100000; ++i)
        {
          double abs_r0_div_abs_r1 = dist_less_than_one(engine);
          // Only r₀ is real, assume r₀ = 1.
          double abs_r1 = 1.0 / abs_r0_div_abs_r1;
          std::uniform_real_distribution<double> real_dist(-abs_r1, abs_r1);
          double alpha = real_dist(engine);
          double beta = std::sqrt(abs_r1 * abs_r1 - alpha * alpha);
          roots(1.0, {alpha, -beta}, {alpha, beta}, graph, point_styles);
        }
        break;
      case 2:   // c
        // r₂ is real, r₀ = α - 𝕚 β, r₁ = α + 𝕚 β.
        for (int i = 0; i < 100000; ++i)
        {
          double abs_r1_div_abs_r2 = dist_less_than_one(engine);
          // Only r₂ is real, assume r₂ = 1.
          double abs_r1 = abs_r1_div_abs_r2;
          std::uniform_real_distribution<double> real_dist(-abs_r1, abs_r1);
          double alpha = real_dist(engine);
          double beta = std::sqrt(abs_r1 * abs_r1 - alpha * alpha);
          roots({alpha, -beta}, {alpha, beta}, 1.0, graph, point_styles);
        }
        break;
#endif
    }
  }

  // Wait with destruction of the window.
  graph.wait_for_keypress();
}
