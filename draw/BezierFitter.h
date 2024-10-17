#pragma once

#define CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS 1

#include "MultiRegion.h"
#include "Line.h"
#include "cairowindow/BezierCurve.h"
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
#include "Point.h"
#endif

namespace cairowindow::draw {

class BezierFitter : public MultiRegion
{
 private:
  std::vector<plot::BezierCurve> plot_bezier_curves_;
  LineStyle line_style_;
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
  PointStyle point_style_;
  std::vector<plot::Point> plot_bezier_curve_points_;   // Begin points of each BezierCurve in plot_bezier_curves_
                                                        // plus the End point of the last BezierCurve in plot_bezier_curves_.
#endif

 public:
  BezierFitter(std::vector<cairowindow::BezierCurve> const& result, LineStyle const& line_style
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
      , PointStyle const& point_style = PointStyle{}
#endif
      ) :
    MultiRegion(line_style.line_color(), line_style.line_width()),
    line_style_(line_style)
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
    , point_style_(point_style)
#endif
  {
    plot_bezier_curves_.reserve(result.size());
    for (cairowindow::BezierCurve const& bezier_curve : result)
    {
      plot_bezier_curves_.emplace_back(bezier_curve);
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
      plot_bezier_curve_points_.emplace_back(bezier_curve.P(0.0));
#endif
    }
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
    if (!result.empty())
      plot_bezier_curve_points_.emplace_back(result.back().P(1.0));
#endif
  }

  // Accessors.
  std::vector<plot::BezierCurve> const& plot_bezier_curves() const { return plot_bezier_curves_; }
  LineStyle line_style() const { return line_style_; }
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
  std::vector<plot::Point> const& plot_bezier_curve_points() const { return plot_bezier_curve_points_; }
  PointStyle point_style() const { return point_style_; }
#endif

 private:
  void draw_regions_on(Layer* layer) override;
};

} // namespace cairowindow::draw
