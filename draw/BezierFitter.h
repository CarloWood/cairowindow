#pragma once

#include "MultiRegion.h"
#include "Line.h"
#include "cairowindow/BezierCurve.h"

namespace cairowindow::draw {

class BezierFitter : public MultiRegion
{
 private:
  std::vector<plot::BezierCurve> plot_bezier_curves_;
  LineStyle line_style_;

 public:
  BezierFitter(std::vector<cairowindow::BezierCurve> const& result, LineStyle const& line_style) :
    MultiRegion(line_style.line_color(), line_style.line_width()), line_style_(line_style)
  {
    plot_bezier_curves_.reserve(result.size());
    for (cairowindow::BezierCurve const& bezier_curve : result)
      plot_bezier_curves_.emplace_back(bezier_curve);
  }

  // Accessors.
  std::vector<plot::BezierCurve> const& plot_bezier_curves() const { return plot_bezier_curves_; }
  LineStyle line_style() const { return line_style_; }

 private:
  void draw_regions_on(Layer* layer) override;
};

} // namespace cairowindow::draw
