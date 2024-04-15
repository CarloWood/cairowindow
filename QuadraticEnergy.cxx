#include "sys.h"
#include "QuadraticEnergy.h"
#include "symbolic/symbolic.h"

namespace cairowindow::autodiff {

// After calling BezierCurve::quadratic_from(double v0qa, double v1qa) (and that returned true),
// the following function can be called to (redo) the calculations for the arc length;
double QuadraticEnergy::arc_length(double v0qa, double v1qa)
{
  v0qa_ = v0qa;
  v1qa_ = v1qa;
  return arc_length_.evaluate();
}

double QuadraticEnergy::stretching_energy(double v0qa, double v1qa)
{
  v0qa_ = v0qa;
  v1qa_ = v1qa;
  return stretching_energy_.evaluate();
}

void QuadraticEnergy::print_derivative() const
{
  auto& d = stretching_energy_.derivative(v0qa_);
  Dout(dc::notice, "derivative ∂stretching_energy_/∂v0qa = " << d << " [" << d.definition() << "] = " << d.evaluate());
}

} // namespace cairowindow::autodiff
