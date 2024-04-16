#include "sys.h"
#include "QuadraticEnergy.h"
#include "symbolic/symbolic.h"

namespace cairowindow::autodiff {

void QuadraticEnergy::init_angles(double v0qa, double v1qa)
{
  v0_div_q1_.reset_evaluation();
  v0x_.reset_evaluation();
  v0y_.reset_evaluation();
  v02_.reset_evaluation();
  A0x_.reset_evaluation();
  A0y_.reset_evaluation();
  a02_.reset_evaluation();
  v0_.reset_evaluation();
  a0_.reset_evaluation();
  z_.reset_evaluation();
  s_.reset_evaluation();
  a03_.reset_evaluation();
  za0pa03s_.reset_evaluation();
  za0v0_.reset_evaluation();
  v02a02mz2_.reset_evaluation();
  zpa02pa0s_.reset_evaluation();
  zpv0a0_.reset_evaluation();
  the_log_.reset_evaluation();
  the_enumerator_.reset_evaluation();
  arc_length_.reset_evaluation();
  stretching_energy_.reset_evaluation();
  dz_.reset_evaluation();
  cs_.reset_evaluation();
  a04_.reset_evaluation();
  dz2_.reset_evaluation();
  e_.reset_evaluation();
  e2_.reset_evaluation();
  se_.reset_evaluation();
  i1_.reset_evaluation();
  i0_.reset_evaluation();
  bending_energy_.reset_evaluation();

  symbolic::Symbol const* v0qa_symbol = dynamic_cast<symbolic::Symbol const*>(&v0qa_);
  symbolic::Symbol const* v1qa_symbol = dynamic_cast<symbolic::Symbol const*>(&v1qa_);

  // There is no need to call init_angles when neither angle is a symbol.
  ASSERT(v0qa_symbol || v1qa_symbol);

  if (v0qa_symbol)
    *v0qa_symbol = v0qa;

  if (v1qa_symbol)
    *v1qa_symbol = v1qa;
}

// After calling BezierCurve::quadratic_from(double v0qa, double v1qa) (and that returned true),
// the following function can be called to (redo) the calculations for the arc length;
double QuadraticEnergy::eval_arc_length()
{
  return arc_length_.evaluate();
}

double QuadraticEnergy::eval_stretching_energy()
{
  return stretching_energy_.evaluate();
}

double QuadraticEnergy::eval_bending_energy()
{
  return bending_energy_.evaluate();
}

} // namespace cairowindow::autodiff
