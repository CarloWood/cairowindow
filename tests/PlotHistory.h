#pragma once

#include "History.h"
#include "PlotSample.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Layer.h"

class PlotHistory : public gradient_descent::History
{
 private:
  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> const& layer_;
  cairowindow::draw::PointStyle const& point_style_;
  cairowindow::draw::TextStyle const& label_style_;
  utils::Array<PlotSample, size, gradient_descent::HistoryIndex> plot_samples_;

 public:
  PlotHistory(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer,
      cairowindow::draw::PointStyle const& point_style, cairowindow::draw::TextStyle const& label_style) :
    plot_(plot), layer_(layer), point_style_(point_style), label_style_(label_style) { }

  void add(double w, double Lw, double dLdw, gradient_descent::Scale const& scale, bool& current_is_replacement)
  {
    gradient_descent::HistoryIndex index = gradient_descent::History::add(w, Lw, dLdw, scale, current_is_replacement);
    plot_samples_[index].initialize(&current(),
      plot_.create_point(layer_, point_style_, {w, Lw}),
      plot_.create_text(layer_, label_style_({.position = cairowindow::draw::centered_below}),
            cairowindow::Point{w, Lw}, std::to_string(total_number_of_samples() - 1)));
  }

  cairowindow::plot::Plot& plot() const { return plot_; }
  boost::intrusive_ptr<cairowindow::Layer> const& layer() const { return layer_; }
  cairowindow::draw::PointStyle const& point_style() const { return point_style_; }
  cairowindow::draw::TextStyle const& label_style() const { return label_style_; }
};
