#pragma once

#include "Scale.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Connector.h"
#include "cairowindow/draw/Connector.h"

class PlotParabolaScale
{
 public:
  static cairowindow::draw::ConnectorStyle const s_indicator_style;

 private:
  gradient_descent::Scale const& master_;

  // Used to visualize the Scale:
  cairowindow::plot::Plot* plot_ = nullptr;
  boost::intrusive_ptr<cairowindow::Layer> layer_;
  cairowindow::plot::Connector plot_indicator_;
  cairowindow::plot::Text plot_scale_text_;
  cairowindow::plot::Line plot_vertical_line_through_w_;
  cairowindow::plot::Line plot_vertical_line_through_v_;

  // Temporary curves (used while developing this class).
  cairowindow::plot::BezierFitter plot_old_parabola_;

 public:
  PlotParabolaScale(gradient_descent::Scale const& master, cairowindow::plot::Plot& plot,
      boost::intrusive_ptr<cairowindow::Layer> const& layer) :
    master_(master), plot_(&plot), layer_(layer) { }

  void draw_indicators()
  {
    ASSERT(master_.has_sample());

    double x1 = master_.parabola().vertex_x();
    double x2 = master_.edge_sample_w();
    double scale_y = plot_->yrange().min() + 0.5 * plot_->yrange().size();
    plot_indicator_ = cairowindow::plot::Connector{{x1, scale_y}, {x2, scale_y},
        cairowindow::Connector::open_arrow, cairowindow::Connector::open_arrow};
    plot_->add_connector(layer_, s_indicator_style, plot_indicator_);
    plot_scale_text_ = plot_->create_text(layer_, {{.position = cairowindow::draw::centered_above, .offset = 2.0}},
        cairowindow::Point{(x1 + x2) / 2, scale_y}, "scale");

    plot_vertical_line_through_v_ = cairowindow::plot::Line{{master_.parabola().vertex_x(), scale_y}, cairowindow::Direction::up};
    plot_->add_line(layer_, s_indicator_style({.dashes = {3.0, 3.0}}), plot_vertical_line_through_v_);

    plot_vertical_line_through_w_ = cairowindow::plot::Line{{master_.edge_sample_w(), scale_y}, cairowindow::Direction::up};
    plot_->add_line(layer_, s_indicator_style, plot_vertical_line_through_w_);
  }

  void erase_indicators()
  {
    // Erase old drawings.
    plot_indicator_.reset();
    plot_scale_text_.reset();
    plot_vertical_line_through_w_.reset();
    plot_vertical_line_through_v_.reset();
    plot_old_parabola_.reset();
  }

  // For use with draw_old_parabola.
  math::QuadraticPolynomial get_parabola() const
  {
    return master_.parabola();
  }

  gradient_descent::Scale const& scale() const
  {
    return master_;
  }

  // Draws and caches a plot of the old parabola; which has to be passed as an argument.
  void draw_old_parabola(math::QuadraticPolynomial const& old_parabola)
  {
    // Draw the old parabola, for debugging purposes.
    using namespace cairowindow;
    plot_old_parabola_.solve([&](double w) -> Point { return {w, old_parabola(w)}; }, plot_->viewport());
    plot_->add_bezier_fitter(layer_, {{.line_color = color::light_red, .line_width = 1.0}}, plot_old_parabola_);
  }
};
