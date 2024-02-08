#pragma once

#include "Point.h"
#include "Range.h"
#include "BezierCurve.h"
#include "IntersectRectangle.h"
#include <functional>
#include <vector>
#include <memory>

namespace cairowindow {

class BezierFitter
{
 private:
  std::vector<BezierCurve> result_;
#ifdef CWDEBUG
  int depth_;
#endif

 public:
  BezierFitter() = default;

  // func     : the parametric function that must be fitted: takes the parameter (t) and returns a Point.
  // domain   : the minimum and maximum values that will be passed to func_ (the domain of t).
  // viewport : Bézier segments that fall entirely outside of this viewport will be discarded.
  // tolerance: the smallest deviation from the true function value allowed in the Bézier curve output.
  void solve(std::function<Point(double)>&& func, Range const& domain, Rectangle const& viewport, double tolerance);

  // Get a pointer to the result vector.
  std::vector<BezierCurve> const& result() const { return result_; }

 private:
  void solve(std::function<Point(double)> const& func, IntersectRectangle const& viewport,
      double tolerance, double t0, double t6, Vector P0, Vector P3, Vector P6);
};

namespace draw {
class BezierFitter;
} // namespace draw

namespace plot {
class Plot;

class BezierFitter : public cairowindow::BezierFitter
{
 public:
  using cairowindow::BezierFitter::BezierFitter;
  BezierFitter(cairowindow::BezierFitter&& bezier_fitter) : cairowindow::BezierFitter(std::move(bezier_fitter)) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::BezierFitter> draw_object_;
};

} // namespace plot
} // namespace cairowindow
