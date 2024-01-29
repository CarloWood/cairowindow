#pragma once

#include "Point.h"
#include "Range.h"
#include "IntersectRectangle.h"
#include "Vector.h"
#include "utils/has_print_on.h"
#include <functional>
#include <vector>

namespace cairowindow {
using utils::has_print_on::operator<<;

class BezierCurve
{
 private:
  Vector P0_;
  Vector P1_;
  Vector P2_;
  Vector P3_;

 public:
  BezierCurve(Vector P0, Vector P1, Vector P2, Vector P3) : P0_(P0), P1_(P1), P2_(P2), P3_(P3) { }

  Vector P0() const { return P0_; }
  Vector P1() const { return P1_; }
  Vector P2() const { return P2_; }
  Vector P3() const { return P3_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

class BezierFitter
{
 private:
  std::function<cairowindow::Point(double)> func_;      // The parametric function that must be fitted: takes the parameter (t) and returns a Point.
  Range domain_;                                        // The minimum and maximum values that will be passed to func_ (the domain of t).
  cairowindow::IntersectRectangle viewport_;            // Bézier segments that fall entirely outside of this viewport will be discarded.
  double tolerance_;                                    // The smallest deviation from the true function value allowed in the Bézier curve output.
#ifdef CWDEBUG
  mutable int depth_;
#endif

 public:
  BezierFitter(std::function<cairowindow::Point(double)>&& func, Range const& domain, cairowindow::Rectangle const& viewport, double tolerance) :
    func_(func), domain_(domain), viewport_(viewport), tolerance_(tolerance) { }

  std::vector<BezierCurve> solve() const;

 private:
  void solve(std::vector<BezierCurve>& out, double t0, double t6, Vector P0, Vector P3, Vector P6) const;
};

} // namespace cairowindow
