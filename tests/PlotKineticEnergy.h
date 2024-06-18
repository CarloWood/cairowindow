#pragma once

#include "KineticEnergy.h"
#include "cairowindow/draw/Connector.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Line.h"
#include "cairowindow/Text.h"

class PlotKineticEnergy : public gradient_descent::KineticEnergy
{
  static cairowindow::draw::ConnectorStyle const s_indicator_style;

 private:
  // Used to visualize the KineticEnergy:
  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> layer_;
  cairowindow::plot::Line plot_horizontal_line_;
  cairowindow::plot::Text plot_energy_text_;

 private:
  void draw_indicators()
  {
    plot_horizontal_line_ = cairowindow::plot::Line{{0.0, max_Lw_}, cairowindow::Direction::right};
    plot_.add_line(layer_, s_indicator_style, plot_horizontal_line_);

    plot_energy_text_ = plot_.create_text(layer_, {{.position = cairowindow::draw::centered_above, .offset = 2.0}},
        cairowindow::Point{0.5 * (plot_.xrange().min() + plot_.xrange().max()), max_Lw_}, "energy");
  }

 public:
  PlotKineticEnergy(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer, double Lw) :
    gradient_descent::KineticEnergy(Lw), plot_(plot), layer_(layer) { }

  void set(double max_Lw, double Lw)
  {
    gradient_descent::KineticEnergy::set(max_Lw, Lw);
    draw_indicators();
  }

  bool maybe_update(double new_Lw)
  {
    if (!gradient_descent::KineticEnergy::maybe_update(new_Lw))
      return false;
    draw_indicators();
    return true;
  }

  void update(double new_Lw)
  {
    gradient_descent::KineticEnergy::update(new_Lw);
    draw_indicators();
  }
};

