#include "sys.h"
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <complex>
#include <cassert>
#include <random>
#include <array>
#include "cairowindow/QuickGraph.h"
#include "debug.h"

void roots(std::complex<double> r0, std::complex<double> r1, std::complex<double> r2, cairowindow::QuickGraph& graph, std::array<cairowindow::draw::PointStyle, 35> const& point_styles)
{
  double c2 = std::real(-(r0 + r1 + r2));
  double c1 = std::real(r0 * r1 + r0 * r2 + r1 * r2);
  double c0 = std::real(-r0 * r1 * r2);

  double f2m30 = std::sqrt(std::abs(c0 * c0 / (c1 * c1 * c1)));

  double y = std::log10(std::abs(r0) / std::abs(r1));

  if (!std::isfinite(y) || y < -4.0)
    return;

  //      |
  //  Q0  |  Q1
  //------+------ <-- y=-2
  //  Q2  |  Q3
  //      |
  //      ^
  //      |
  //     x=-2

  int quadrant = (f2m30 > 0.01 ? 1 : 0) | (y < -2 ? 2 : 0);

#if 0
  if (quadrant == 1 || quadrant == 3)
  {
    int count = 0;
    for (int e0 = -3; e0 <= 3; ++e0)
      for (int e1 = -3; e1 <= 3; ++e1)
      {
        int e2 = -(3 * e0 + 2 * e1);
        if (e2 < -3 || e2 > 3)
          continue;
//        if (e0 == 2 && e1 == -3)
//          continue;
        double f = std::pow(c0, e0) * std::pow(c1, e1) * std::pow(c2, e2);
        double x = 0.5 * std::log10(f);

        if (std::isfinite(x) && -8.0 < x && x < 4.0)
          graph.add_point({x, -4.0 * count / 17.0 - quadrant / 40.0}, point_styles[count + 9 * (quadrant - 1)]);

        ++count;
      }
  }
#else
  double x = std::log10(f2m30);
  if (std::isfinite(x) && -8.0 < x && x < 4.0)
    graph.add_point({x, y}, point_styles[0]);
#endif
}

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

  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  std::uniform_real_distribution<double> dist_less_than_one(0.0001, 1.0);

  for (;;)
  for (int case_abc = 1; case_abc < 3; ++case_abc)
  {
    switch (case_abc)
    {
      case 0:   // a
      {
        for (int i = 0; i < 100000; ++i)
        {
          // All three roots are reals.
          double r1_div_r2 = dist(engine);
          double r0_div_r1 = dist(engine);
          if (std::abs(r1_div_r2) < 0.0001)
            continue;
          roots(r0_div_r1, 1.0, 1.0 / r1_div_r2, graph, point_styles);
        }
        break;
      }
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
    }
  }

  // Wait with destruction of the window.
  graph.wait_for_keypress();
}
