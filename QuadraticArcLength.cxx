#include "sys.h"
#include "QuadraticArcLength.h"
#include "symbolic/symbolic.h"

namespace cairowindow::autodiff {

// After calling BezierCurve::quadratic_from(double v0qa, double v1qa) (and that returned true),
// the following function can be called to (redo) the calculations for the arc length;
double QuadraticArcLength::quadratic_arc_length(double v0qa, double v1qa)
{
  v0qa_ = v0qa;
  v1qa_ = v1qa;
  return evaluate(quadratic_arc_length_);
}

} // namespace cairowindow::autodiff
