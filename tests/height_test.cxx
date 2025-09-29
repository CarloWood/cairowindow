#include "sys.h"
#include "gradient_descent2/ExtremeType.h"
#include "gradient_descent2/AnalyzedCubic.h"
#include "math/CubicPolynomial.h"
#include <array>
#include <iomanip>
#include <limits>
#include "debug.h"

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  using namespace gradient_descent;

  // Generate a cubic that has its extrema at -3 and 14.
  constexpr int si = -3;
  constexpr int ti = 14;
  int const ai = 42;

  // By default assume si is the maximum.
  // Lets have two points left of the maximum, three in between, and two on the right of the minimum.
  static_assert(ti > si, "ti must be larger than si");
  static_assert((si + ti) / 2 > si + 1, "ti - si is too small");
  static_assert((si + ti) / 2 < ti - 2, "ti - si is too small");
  std::array<int, 9> xs = { si - 2, si - 1, si, si + 1, (si + ti) / 2, ti - 2, ti, ti + 2, ti + 3 };

  constexpr int i_max = 2;    //    ^^
  constexpr int i_min = 6;    //                                               ^^

  // The inflection point is at
  constexpr double inflection_point = 0.5 * (si + ti);

  for (int maximum_last = 0; maximum_last <= 1; ++ maximum_last)
  {
    bool inverted = maximum_last;                           // True if ti is the maximum.
    int di = inverted ? -2 : 2;
    int bi = 3 * di * si * ti;
    int two_ci = -3 * di * (si + ti);

    double a = ai;
    double b = bi;
    double c = 0.5 * two_ci;
    double d = di;

    math::CubicPolynomial<double> g(a, b, c, d);

    Dout(dc::notice, "g(x) = " << g);

    std::array<double, 2> extrema;
    Debug(dc::notice.off());
    int number_of_extrema = g.get_extrema(extrema, false);
    Debug(dc::notice.on());
    ASSERT(number_of_extrema == 2);
    ASSERT(g(extrema[0]) < g(extrema[1]));

    std::array<ExtremeType, 2> extreme_type = { ExtremeType::minimum, ExtremeType::maximum };

    for (int i = 0; i <= 1; ++i)
    {
      Dout(dc::notice, "  " << extreme_type[i] << " at " << extrema[i] << ": " << g(extrema[i]));
      AnalyzedCubic acubic;
      acubic.initialize(g, extreme_type[i]);

      double extreme_w = acubic.get_extreme();
      Dout(dc::notice, "    extreme at: " << extreme_w);

      constexpr double step = 0.001;
      for (double dw = -step; dw < 1.5 * step; dw += step)
      {
        Dout(dc::notice, "    height at extreme + " << dw << ": " << std::setprecision(std::numeric_limits<double>::digits10) <<
            acubic.height(extreme_w + dw, g[3]));
        Dout(dc::notice, "    actual height: " << std::setprecision(std::numeric_limits<double>::digits10) <<
            std::abs(g(extreme_w + dw) - g(extreme_w)));
      }
    }
  }
}
