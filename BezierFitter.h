#pragma once

#include "BezierCurve.h"
#include "Point.h"
#include "Vector.h"
#include "Rectangle.h"
#include "Range.h"
#include "IntersectRectangle.h"
#include <functional>
#include <vector>
#include <memory>

namespace cairowindow {

// Whether the distance of a point to a fitted Bezier curve is measured horizontally, vertically or perpendicular to the curve.
enum class Orientation
{
  horizontal,
  vertical,     // This is what one should typically use for normal functions.
  perpendicular // This should only be used when rotations of the curve make sense (x and y are interchangable and the draw ratio is 1:1).
};

#ifdef CWDEBUG
std::string to_string(Orientation orientation);
inline std::ostream& operator<<(std::ostream& os, Orientation orientation)
{
  return os << to_string(orientation);
}
#endif

class BezierFitter
{
 private:
  std::vector<BezierCurve> result_;
  int depth_;

 public:
  BezierFitter() = default;

  void solve(std::function<void(Point p, Vector v)> const& draw_line,
      std::function<Point(double)>&& P, std::function<Vector(double)>&& V,
      Range const& domain, Rectangle const& viewport, double fraction, Orientation orientation);

  // func     : the parametric function that must be fitted: takes the parameter (t) and returns a Point.
  // domain   : the minimum and maximum values that will be passed to func_ (the domain of t).
  // viewport : Bézier segments that fall entirely outside of this viewport will be discarded.
  //            Normally you'd pass plot.viewport() here.
  // tolerance: the smallest deviation from the true function value allowed in the Bézier curve output;
  //            when -1.0 is passed the default 1e-5 * viewport.height() will be used.
  void solve(std::function<Point(double)>&& func, Range const& domain, Rectangle const& viewport, double tolerance = -1.0);

  // Convenience constructor that combine the construction with the call to solve.
  BezierFitter(std::function<Point(double)>&& func, Range const& domain, Rectangle const& viewport, double tolerance = -1.0)
  {
    solve(std::move(func), domain, viewport, tolerance);
  }

  // This signature can be used, for example, when `func` returns a Point that simply passes its argument back as the x-coordinate of the Point.
  void solve(std::function<Point(double)>&& func, Rectangle const& viewport, double tolerance = -1.0)
  {
    solve(std::move(func), {viewport.offset_x(), viewport.offset_x() + viewport.width()}, viewport, tolerance);
  }

  // Same
  BezierFitter(std::function<Point(double)>&& func, Rectangle const& viewport, double tolerance = -1.0)
  {
    solve(std::move(func), viewport, tolerance);
  }

  // Get a pointer to the result vector.
  std::vector<BezierCurve> const& result() const { return result_; }

 private:
  void solve(std::function<void(Point p, Vector v)> const& draw_line,
      std::function<Point(double)> const& P, std::function<Vector(double)> const& V,
      IntersectRectangle<csid::plot> const& viewport, double fraction, Orientation orientation,
      double t0, double t1, BezierCurve::Point P0, BezierCurve::Vector T0, BezierCurve::Point Pg, BezierCurve::Point P1, BezierCurve::Vector T1);

  void solve(std::function<Point(double)> const& func, IntersectRectangle<csid::plot> const& viewport,
      double tolerance, double t0, double t6, BezierCurve::Vector P0, BezierCurve::Vector P3, BezierCurve::Vector P6);
};

namespace draw {
class BezierFitter;
} // namespace draw

namespace plot {
class Plot;

class BezierFitter : public cairowindow::BezierFitter
{
 public:
  explicit BezierFitter(cairowindow::BezierFitter&& bezier_fitter) : cairowindow::BezierFitter(std::move(bezier_fitter)) { }
  using cairowindow::BezierFitter::BezierFitter;

  void reset()
  {
    draw_object_.reset();
  }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::BezierFitter> draw_object_;
};

} // namespace plot
} // namespace cairowindow
