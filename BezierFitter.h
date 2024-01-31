#pragma once

#include "Point.h"
#include "Range.h"
#include "BezierCurve.h"
#include "IntersectRectangle.h"
#include <functional>
#include <vector>

namespace cairowindow {

class BezierFitter
{
 private:
  std::function<Point(double)> func_;           // The parametric function that must be fitted: takes the parameter (t) and returns a Point.
  Range domain_;                                // The minimum and maximum values that will be passed to func_ (the domain of t).
  IntersectRectangle viewport_;                 // Bézier segments that fall entirely outside of this viewport will be discarded.
  double tolerance_;                            // The smallest deviation from the true function value allowed in the Bézier curve output.
#ifdef CWDEBUG
  mutable int depth_;
#endif

 public:
  BezierFitter(std::function<Point(double)>&& func, Range const& domain, Rectangle const& viewport, double tolerance) :
    func_(func), domain_(domain), viewport_(viewport), tolerance_(tolerance) { }

  std::vector<BezierCurve> solve() const;

 private:
  void solve(std::vector<BezierCurve>& out, double t0, double t6, Vector P0, Vector P3, Vector P6) const;
};

} // namespace cairowindow
