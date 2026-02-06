#pragma once

#include "gradient_descent2/Sample.h"
#include "gradient_descent2/AlgorithmEventType.h"
#include "cairowindow/Point.h"
#include "cairowindow/Text.h"
#include "cairowindow/Connector.h"
#include "cairowindow/BezierFitter.h"
#include "cairowindow/Line.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Direction.h"
#include "cairowindow/draw/Point.h"
#include "cairowindow/draw/Text.h"
#include "cairowindow/draw/Connector.h"
#include "utils/Array.h"

#ifdef CWDEBUG
// A Sample including plot objects: a point and a label.
class PlotSample
{
 private:
  cairowindow::plot::Point P_;
  cairowindow::plot::Text P_label_;

 public:
  PlotSample() = default;
  PlotSample(cairowindow::plot::Point const& point, cairowindow::plot::Text const& label) : P_(point), P_label_(label) { }

  void initialize(cairowindow::plot::Point const& point, cairowindow::plot::Text const& label)
  {
    P_ = point;
    P_label_ = label;
  }

  void replace(cairowindow::plot::Point const& point)
  {
    P_ = point;
  }

  cairowindow::Point const& P() const { return P_; }
  cairowindow::Text const& label() const { return P_label_; }

  std::string debug_label() const { return P_label_.text(); }

  void print_on(std::ostream& os) const
  {
    os << debug_label() << " (at w = " << P_.x() << ")";
  }
};

class AlgorithmEvent
{
 public:
  static constexpr cairowindow::draw::PointStyle point_style_{{.color_index = 31, .filled_shape = 1}};
  static cairowindow::draw::TextStyle const s_label_style;
  static constexpr cairowindow::draw::LineStyle curve_line_style_{{.line_width = 1.0}};
  static cairowindow::draw::ConnectorStyle const s_difference_expected_style;
  static cairowindow::draw::ConnectorStyle const s_indicator_style;

 private:
  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> const& layer_;
  cairowindow::plot::BezierFitter plot_quadratic_approximation_curve_;
  cairowindow::plot::BezierFitter plot_fourth_degree_approximation_curve_;
  cairowindow::plot::BezierFitter plot_quotient_curve_;
  cairowindow::plot::BezierFitter plot_derivative_curve_;
  cairowindow::plot::Connector plot_difference_expected_;
  cairowindow::plot::Line plot_horizontal_line_;
  cairowindow::plot::Text plot_energy_text_;
  std::array<cairowindow::plot::Connector, 3> plot_scale_indicator_;    // 0: left, 1: right, 2: scale
  std::array<cairowindow::plot::Text, 3> plot_scale_text_;              // idem
  std::array<cairowindow::plot::Line, 4> plot_vertical_line_through_w_; // 0: left, 1: right, 2: cp, 3: I
  std::array<cairowindow::plot::BezierFitter, 2> plot_cubic_curve_;     // 0: left, 1: right
  cairowindow::plot::BezierFitter plot_old_cubic_;
  std::array<PlotSample, 32> plot_samples_;                             // Circular buffer for the last 32 added samples.
  int current_{};                                                       // Index into plot_samples_.
  std::vector<PlotSample> plot_local_extremes_;
  cairowindow::plot::Connector plot_current_hdirection_;
  cairowindow::plot::Connector plot_current_left_of_direction_;
  cairowindow::plot::Connector plot_current_right_of_direction_;
  cairowindow::plot::Connector plot_next_jump_point_;

 public:
  AlgorithmEvent(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer) :
    plot_(plot), layer_(layer) { }

  void callback(gradient_descent::AlgorithmEventType const& event);
};

#endif
