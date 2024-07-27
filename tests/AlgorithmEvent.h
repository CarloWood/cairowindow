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
  cairowindow::plot::BezierFitter plot_cubic_approximation_curve_;
  cairowindow::plot::BezierFitter plot_fourth_degree_approximation_curve_;
  cairowindow::plot::BezierFitter plot_quotient_curve_;
  cairowindow::plot::BezierFitter plot_derivative_curve_;
  cairowindow::plot::Connector plot_difference_expected_;
  cairowindow::plot::Line plot_horizontal_line_;
  cairowindow::plot::Text plot_energy_text_;
  std::array<cairowindow::plot::Connector, 3> plot_scale_indicator_;    // 0: left, 1: right, 2: scale
  std::array<cairowindow::plot::Text, 3> plot_scale_text_;              // idem
  std::array<cairowindow::plot::Line, 4> plot_vertical_line_through_w_; // 0: left, 1: right, 2: cp, 3: I
  cairowindow::plot::BezierFitter plot_old_cubic_;
  utils::Array<PlotSample, gradient_descent::History::size, gradient_descent::HistoryIndex> plot_samples_;
  std::vector<PlotSample> plot_local_extremes_;
  cairowindow::plot::Connector plot_current_hdirection_;

 public:
  AlgorithmEvent(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer) :
    plot_(plot), layer_(layer) { }

  void callback(gradient_descent::AlgorithmEventType const& event)
  {
    using namespace cairowindow;
    using namespace gradient_descent;

    if (event.is_a<ResetEventData>())
    {
      plot_approximation_curve_.reset();
      plot_derivative_curve_.reset();
      plot_quotient_curve_.reset();
      plot_fourth_degree_approximation_curve_.reset();
    }
    else if (event.is_a<DifferenceEventData>())
    {
      auto const& data = event.get<DifferenceEventData>();

      // Plot the vertical difference from what we expected to what we got.
      plot_difference_expected_ = cairowindow::plot::Connector{{data.w(), data.expected_Lw()}, {data.w(), data.Lw()},
          cairowindow::Connector::no_arrow, cairowindow::Connector::open_arrow};
      plot_.add_connector(layer_, s_difference_expected_style, plot_difference_expected_);
    }
    else if (event.is_a<FourthDegreeApproximationEventData>())
    {
      auto const& data = event.get<FourthDegreeApproximationEventData>();

      math::Polynomial const& fourth_degree_approximation = data.polynomial();
      plot_fourth_degree_approximation_curve_.solve(
          [&fourth_degree_approximation](double w) -> Point { return {w, fourth_degree_approximation(w)}; }, plot_.viewport());
      plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::teal}), plot_fourth_degree_approximation_curve_);
    }
    else if (event.is_a<DerivativeEventData>())
    {
      auto const& data = event.get<DerivativeEventData>();

      math::Polynomial const& derivative = data.polynomial();
      plot_derivative_curve_.solve([&derivative](double w) -> Point { return {w, derivative(w)}; }, plot_.viewport());
      plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::magenta}), plot_derivative_curve_);
    }
    else if (event.is_a<QuotientEventData>())
    {
      auto const& data = event.get<QuotientEventData>();

      math::Polynomial const& quotient = data.polynomial();
      plot_quotient_curve_.solve([&quotient](double w) -> Point { return {w, 10.0 * quotient(w)}; }, plot_.viewport());
      plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::blue}), plot_quotient_curve_);
    }
    else if (event.is_a<QuadraticPolynomialEventData>())
    {
      auto const& data = event.get<QuadraticPolynomialEventData>();

      math::QuadraticPolynomial const& approximation = data.quadratic_polynomial();
      plot_approximation_curve_.solve([&approximation](double w) -> Point { return {w, approximation(w)}; }, plot_.viewport());
      plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::red}), plot_approximation_curve_);
    }
    else if (event.is_a<KineticEnergyEventData>())
    {
      auto const& data = event.get<KineticEnergyEventData>();

      plot_horizontal_line_ = plot::Line{{0.0, data.max_Lw()}, Direction::right};
      plot_.add_line(layer_, s_indicator_style, plot_horizontal_line_);

      plot_energy_text_ = plot_.create_text(layer_, {{.position = draw::centered_above, .offset = 2.0}},
          Point{0.5 * (plot_.xrange().min() + plot_.xrange().max()), data.max_Lw()}, "energy");
    }
    else if (event.is_a<ScaleDrawEventData>())
    {
      auto const& data = event.get<ScaleDrawEventData>();

      using namespace gradient_descent;
      if (data.result() == ScaleUpdate::initialized ||
          data.result() == ScaleUpdate::towards_cp ||
          data.result() == ScaleUpdate::away_from_cp)
      {
        Scale const& scale = data.scale();
        double const I = scale.inflection_point_w();
        double const cp = scale.critical_point_w();
        double const l = scale.left_edge_sample_w();
        double const r = scale.right_edge_sample_w();
        double prev_x2;
        double const y = plot_.yrange().min() + 0.75 * plot_.yrange().size();
        double scale_y = y;
        for (int side = 0; side < 2; ++side)    // left and right
        {
          double x2 = side == 0 ? l : r;
#if 0   // Do not draw the left/right connectors for now.
          if (side == 1 && std::copysign(1.0, x2 - cp) == std::copysign(1.0, prev_x2 - cp))
            scale_y -= plot_.convert_vertical_offset_from_pixel(15.0);
          prev_x2 = x2;
          plot_scale_indicator_[side] = plot::Connector{{cp, scale_y}, {x2, scale_y}, Connector::open_arrow, Connector::open_arrow};
          plot_.add_connector(layer_, s_indicator_style, plot_scale_indicator_[side]);
          plot_scale_text_[side] = plot_.create_text(layer_, {{.position = draw::centered_above, .offset = 2.0}},
              Point{(cp + x2) / 2, scale_y}, side == 0 ? "left" : "right");
#endif

          // Only draw the vertical line at critical_point_w and inflection point once.
          if (side == 0)
          {
            plot_vertical_line_through_w_[2] = plot::Line{{cp, 0.0}, Direction::up};
            plot_.add_line(layer_, s_indicator_style({.dashes = {3.0, 3.0}}), plot_vertical_line_through_w_[2]);
            plot_vertical_line_through_w_[3] = plot::Line{{I, 0.0}, Direction::up};
            plot_.add_line(layer_, s_indicator_style({.line_color = color::violet, .dashes = {5.0, 2.0}}), plot_vertical_line_through_w_[3]);
          }

          plot_vertical_line_through_w_[side] = plot::Line{{x2, 0.0}, Direction::up};
          plot_.add_line(layer_, s_indicator_style, plot_vertical_line_through_w_[side]);
        }
        // Draw the "scale" connector.
        scale_y = y + plot_.convert_vertical_offset_from_pixel(15.0);
        double x_center;
        Scale::LCRI_type const LCRI{l, cp, r, I};
        switch (Scale::classify(LCRI))
        {
          case 0:       // Weighted average of {I, r} around cp.
            x_center = cp;
            break;
          case 1:       // Weighted average of {l, r} around I.
            x_center = I;
            break;
          case 2:       // Weighted average of {l, I} around cp.
            x_center = cp;
            break;
          case 3:       // Weighted average of {l, r} around cp.
            x_center = cp;
            break;
          case 4:       // r - l
            x_center = 0.5 * (l + r);
            break;
        }
        double width = scale.value();
        plot_scale_indicator_[2] = plot::Connector{{x_center - 0.5 * width, scale_y}, {x_center + 0.5 * width, scale_y},
          Connector::open_arrow, Connector::open_arrow};
        plot_.add_connector(layer_, s_indicator_style, plot_scale_indicator_[2]);
        plot_scale_text_[2] = plot_.create_text(layer_, {{.position = draw::centered_above, .offset = 2.0}}, Point{x_center, scale_y}, "scale");
      }
      if (data.result() == ScaleUpdate::towards_cp)
      {
        auto const& old_cubic = data.old_cubic();
        // Draw the old cubic.
        plot_old_cubic_.solve([&old_cubic](double w) -> Point { return {w, old_cubic(w)}; }, plot_.viewport());
        plot_.add_bezier_fitter(layer_, {{.line_color = color::lightred, .line_width = 1.0}}, plot_old_cubic_);
      }
    }
    else if (event.is_a<ScaleEraseEventData>())
    {
      for (int side = 0; side < plot_scale_indicator_.size(); ++side)
      {
        plot_scale_indicator_[side].reset();
        plot_scale_text_[side].reset();
      }
      for (int i = 0; i < plot_vertical_line_through_w_.size(); ++i)
        plot_vertical_line_through_w_[i].reset();
      plot_old_cubic_.reset();
    }
    else if (event.is_a<HistoryAddEventData>())
    {
      auto const& data = event.get<HistoryAddEventData>();

      plot_samples_[data.index()].initialize(&data.current(),
        plot_.create_point(layer_, point_style_, {data.current().w(), data.current().Lw()}),
        plot_.create_text(layer_, s_label_style({.position = cairowindow::draw::centered_below}),
              cairowindow::Point{data.current().w(), data.current().Lw()}, data.label()));
    }
    else if (event.is_a<CubicPolynomialEventData>())
    {
      auto const& data = event.get<CubicPolynomialEventData>();

      math::CubicPolynomial const& cubic_approximation = data.cubic_polynomial();
      plot_cubic_approximation_curve_.solve([&cubic_approximation](double w) -> Point { return {w, cubic_approximation(w)}; }, plot_.viewport());
      plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::magenta}), plot_cubic_approximation_curve_);
    }
    else if (event.is_a<NewLocalExtremeEventData>())
    {
      auto const& data = event.get<NewLocalExtremeEventData>();

      plot_local_extremes_.emplace_back(&data.local_extreme(),
          plot_.create_point(layer_, point_style_({.color_index = 2}), {data.local_extreme().w(), data.local_extreme().Lw()}),
          plot_.create_text(layer_, s_label_style({.position = cairowindow::draw::centered_above}),
            cairowindow::Point{data.local_extreme().w(), data.local_extreme().Lw()}, data.label()));
    }
    else if (event.is_a<HDirectionKnownEventData>())
    {
      auto const& data = event.get<HDirectionKnownEventData>();
      double x = data.local_extreme().w();
      double y = data.local_extreme().Lw();

      plot_current_hdirection_ = plot::Connector{{x, y},
        {x + static_cast<int>(data.hdirection()) * plot_.convert_horizontal_offset_from_pixel(25.0), y},
        Connector::no_arrow, Connector::open_arrow};
      plot_.add_connector(layer_, s_indicator_style({.line_color = color::red}), plot_current_hdirection_);
    }
    else
      // Missing implementation.
      ASSERT(false);
  }
};

#endif
