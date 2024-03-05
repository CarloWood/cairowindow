#pragma once

#include "MultiRegion.h"
#include "BezierCurve.h"
#include <memory>
#include <vector>

namespace cairowindow::draw {

class Curve : public MultiRegion
{
 private:
  std::vector<std::shared_ptr<BezierCurve>> bezier_curves_;

  // Implementation of MultiRegion.
  void draw_regions_on(Layer* layer) override;

 public:
  Curve(BezierCurveStyle const& bezier_curve_style) :
    MultiRegion(bezier_curve_style.line_color(), bezier_curve_style.line_width()) { }

  Curve(std::vector<std::shared_ptr<BezierCurve>>&& bezier_curves, BezierCurveStyle const& bezier_curve_style) :
    MultiRegion(bezier_curve_style.line_color(), bezier_curve_style.line_width()), bezier_curves_(std::move(bezier_curves)) { }

  std::vector<std::shared_ptr<BezierCurve>> const& bezier_curves() const { return bezier_curves_; }
  std::vector<std::shared_ptr<BezierCurve>>& bezier_curves() { return bezier_curves_; }
};

} // namespace cairowindow::draw
