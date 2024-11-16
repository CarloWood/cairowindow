#include "sys.h"
#include "math/bracket_zero.h"
#include "debug.h"
#include <cmath>

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  // The mathematical solution of (1e16 + x) * 1.001 - 1e16 = 0 is x = -9990009990009.99000999000999000999000999...
  // However, since the 1.001 is really stored as 1.0009999999999998899, the best mathematical solution is x = -9990009990009.9911087913...
  //
  // The best double representation of both is -9990009990009.990 (Δ = ±0.002)
  //
  // However, the best representation of the intermediate value (1e16 + x) = 9990009990009990 (Δ = ±2)
  // where the only value of y for which y * 1.001 returns 10000000000000000 is 9990009990009992 (Δ = ±2).
  // Hence that math::bracket_zero finds 9990009990009992 for the intermediate value and thus finds x = -9990009990008.
  //
  // Given that the resolution can't be better than Δ = ±2 because of the intermediate value,
  // the best answer would have been -9990009990010 (Δ = ±2)! So why isn't it?
  //
  // The reason is that for x = -9990009990010, 1e16 + x = 9990009990009990 (which is exact) but
  // 1.001 * 9990009990009990 = (4508103226997866 / 4503599627370496) * (4995004995004996 * 2) =
  // 22517998136852477499959524340670 / 2251799813685248 = 9999999999999998.89 = 9999999999999998 (Δ = ±2).
  // and thus f(-9990009990010) returns -2, not zero.
  //
  // That the influence of the precision of 1.001 is profound can be seen by using `std::nextafter(1.001, 2.0)` instead
  // of 1.001 when f(-9990009990010) return +2!
  //
  auto f = [](double x){ return (1e16 + x) * 1.001 - 1e16; };

  double result = math::bracket_zero(-1e20, 1e20, f);

  Dout(dc::notice, std::setprecision(18) << "result = " << result);
}
