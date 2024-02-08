#include "sys.h"
#include "BezierFitter.h"

namespace cairowindow::draw {

void BezierFitter::draw_regions_on(Layer* layer)
{
  // Nothing to do here, because draw::BezierFitter stores a vector of plot::BezierCurve's,
  // which do their own drawing when added to the plot.
}

} // namespace cairowindow::draw
