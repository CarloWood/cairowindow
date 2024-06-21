#pragma once

#include "gradient_descent/Sample.h"
#include "gradient_descent/History.h"
#include "cairowindow/Point.h"
#include "cairowindow/Text.h"
#include "cairowindow/Connector.h"
#include "cairowindow/BezierFitter.h"
#include "cairowindow/Line.h"
#include "cairowindow/Plot.h"
#include "cairowindow/draw/Point.h"
#include "cairowindow/draw/Text.h"
#include "cairowindow/draw/Connector.h"
#include "utils/Array.h"

#ifdef CWDEBUG
// A Sample including plot objects: a point and a label.
class PlotSample
{
 private:
  gradient_descent::Sample const* master_{};

  mutable cairowindow::plot::Point P_;
  mutable cairowindow::plot::Text P_label_;

 public:
  // Required for History.
  PlotSample() = default;

  // Constructor used by LocalExtreme.
  PlotSample(gradient_descent::Sample const* master, cairowindow::plot::Point const& point, cairowindow::plot::Text const& label) :
    master_(master), P_(point), P_label_(label) { }

  // Required for History.
  void initialize(gradient_descent::Sample const* master, cairowindow::plot::Point const& point, cairowindow::plot::Text const& label)
  {
    master_ = master;
    P_ = point;
    P_label_ = label;
  }

  cairowindow::Point const& P() const { return P_; }
  cairowindow::Text const& label() const { return P_label_; }

  gradient_descent::Sample const* sample() const { return master_; }

  double w() const { return master_->w(); }
  double Lw() const { return master_->Lw(); }
  double dLdw() const { return master_->dLdw(); }

  std::string debug_label() const { return P_label_.text(); }

  void print_on(std::ostream& os) const
  {
    os << debug_label() << " (at w = " << master_->w() << ")";
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
  cairowindow::plot::BezierFitter plot_approximation_curve_;
  cairowindow::plot::BezierFitter plot_fourth_degree_approximation_curve_;
  cairowindow::plot::BezierFitter plot_quotient_curve_;
  cairowindow::plot::BezierFitter plot_derivative_curve_;
  cairowindow::plot::Connector plot_difference_expected_;
  cairowindow::plot::Line plot_horizontal_line_;
  cairowindow::plot::Text plot_energy_text_;
  cairowindow::plot::Connector plot_scale_indicator_;
  cairowindow::plot::Text plot_scale_text_;
  cairowindow::plot::Line plot_vertical_line_through_w_;
  cairowindow::plot::Line plot_vertical_line_through_v_;
  cairowindow::plot::BezierFitter plot_old_parabola_;
  utils::Array<PlotSample, gradient_descent::History::size, gradient_descent::HistoryIndex> plot_samples_;

 public:
  AlgorithmEvent(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer) :
    plot_(plot), layer_(layer) { }

  void callback(gradient_descent::AlgorithmEventType const& event)
  {
    using namespace cairowindow;

    if (event.is_a<gradient_descent::ResetEventData>())
    {
      plot_approximation_curve_.reset();
      plot_derivative_curve_.reset();
      plot_quotient_curve_.reset();
      plot_fourth_degree_approximation_curve_.reset();
    }
    else if (event.is_a<gradient_descent::DifferenceEventData>())
    {
      auto const& data = event.get<gradient_descent::DifferenceEventData>();

      // Plot the vertical difference from what we expected to what we got.
      plot_difference_expected_ = cairowindow::plot::Connector{{data.w(), data.expected_Lw()}, {data.w(), data.Lw()},
          cairowindow::Connector::no_arrow, cairowindow::Connector::open_arrow};
      plot_.add_connector(layer_, s_difference_expected_style, plot_difference_expected_);
    }
    else if (event.is_a<gradient_descent::FourthDegreeApproximationEventData>())
    {
      auto const& data = event.get<gradient_descent::FourthDegreeApproximationEventData>();

      math::Polynomial const& fourth_degree_approximation = data.polynomial();
      plot_fourth_degree_approximation_curve_.solve(
          [&fourth_degree_approximation](double w) -> Point { return {w, fourth_degree_approximation(w)}; }, plot_.viewport());
      plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::teal}), plot_fourth_degree_approximation_curve_);
    }
    else if (event.is_a<gradient_descent::DerivativeEventData>())
    {
      auto const& data = event.get<gradient_descent::DerivativeEventData>();

      math::Polynomial const& derivative = data.polynomial();
      plot_derivative_curve_.solve([&derivative](double w) -> Point { return {w, derivative(w)}; }, plot_.viewport());
      plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::magenta}), plot_derivative_curve_);
    }
    else if (event.is_a<gradient_descent::QuotientEventData>())
    {
      auto const& data = event.get<gradient_descent::QuotientEventData>();

      math::Polynomial const& quotient = data.polynomial();
      plot_quotient_curve_.solve([&quotient](double w) -> Point { return {w, 10.0 * quotient(w)}; }, plot_.viewport());
      plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::blue}), plot_quotient_curve_);
    }
    else if (event.is_a<gradient_descent::QuadraticPolynomialEventData>())
    {
      auto const& data = event.get<gradient_descent::QuadraticPolynomialEventData>();

      math::QuadraticPolynomial const& approximation = data.quadratic_polynomial();
      plot_approximation_curve_.solve([&approximation](double w) -> Point { return {w, approximation(w)}; }, plot_.viewport());
      plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::red}), plot_approximation_curve_);
    }
    else if (event.is_a<gradient_descent::KineticEnergyEventData>())
    {
      auto const& data = event.get<gradient_descent::KineticEnergyEventData>();

      plot_horizontal_line_ = plot::Line{{0.0, data.max_Lw()}, Direction::right};
      plot_.add_line(layer_, s_indicator_style, plot_horizontal_line_);

      plot_energy_text_ = plot_.create_text(layer_, {{.position = draw::centered_above, .offset = 2.0}},
          Point{0.5 * (plot_.xrange().min() + plot_.xrange().max()), data.max_Lw()}, "energy");
    }
    else if (event.is_a<gradient_descent::ScaleDrawEventData>())
    {
      auto const& data = event.get<gradient_descent::ScaleDrawEventData>();

      using namespace gradient_descent;
      if (data.result() == ScaleUpdate::initialized ||
          data.result() == ScaleUpdate::towards_vertex ||
          data.result() == ScaleUpdate::away_from_vertex)
      {
        double x1 = data.x1();
        double x2 = data.x2();
        double scale_y = plot_.yrange().min() + 0.5 * plot_.yrange().size();
        plot_scale_indicator_ = plot::Connector{{x1, scale_y}, {x2, scale_y},
            Connector::open_arrow, Connector::open_arrow};
        plot_.add_connector(layer_, s_indicator_style, plot_scale_indicator_);
        plot_scale_text_ = plot_.create_text(layer_, {{.position = draw::centered_above, .offset = 2.0}},
            Point{(x1 + x2) / 2, scale_y}, "scale");

        plot_vertical_line_through_v_ = plot::Line{{x1, scale_y}, Direction::up};
        plot_.add_line(layer_, s_indicator_style({.dashes = {3.0, 3.0}}), plot_vertical_line_through_v_);

        plot_vertical_line_through_w_ = plot::Line{{x2, scale_y}, Direction::up};
        plot_.add_line(layer_, s_indicator_style, plot_vertical_line_through_w_);
      }
      if (data.result() == ScaleUpdate::towards_vertex)
      {
        auto const& old_parabola = data.old_parabola();
        // Draw the old parabola.
        plot_old_parabola_.solve([&old_parabola](double w) -> Point { return {w, old_parabola(w)}; }, plot_.viewport());
        plot_.add_bezier_fitter(layer_, {{.line_color = color::light_red, .line_width = 1.0}}, plot_old_parabola_);
      }
    }
    else if (event.is_a<gradient_descent::ScaleEraseEventData>())
    {
      plot_scale_indicator_.reset();
      plot_scale_text_.reset();
      plot_vertical_line_through_w_.reset();
      plot_vertical_line_through_v_.reset();
      plot_old_parabola_.reset();
    }
    else if (event.is_a<gradient_descent::HistoryAddEventData>())
    {
      auto const& data = event.get<gradient_descent::HistoryAddEventData>();

      plot_samples_[data.index()].initialize(&data.current(),
        plot_.create_point(layer_, point_style_, {data.current().w(), data.current().Lw()}),
        plot_.create_text(layer_, s_label_style({.position = cairowindow::draw::centered_below}),
              cairowindow::Point{data.current().w(), data.current().Lw()}, data.label()));
    }
  }
};

#endif
