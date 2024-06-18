#pragma once

#include "LocalExtreme.h"
#include "Approximation.h"
#include "PlotParabolaScale.h"
#include "PlotSample.h"
#include "cairowindow/Point.h"

struct PlotLocalExtreme : gradient_descent::LocalExtreme
{
  PlotParabolaScale plot_approximation_parabola_scale_;         // Sibling of gradient_descent::LocalExtreme::approximation_
  PlotSample plot_vertex_sample_;                               // Sibling of gradient_descent::LocalExtreme::vertex_sample_

  PlotLocalExtreme(PlotHistory const& history, gradient_descent::Approximation const& approximation, double energy) :
    gradient_descent::LocalExtreme(history.current(), approximation, energy),
    plot_approximation_parabola_scale_(approximation_.parabola_scale(), history.plot(), history.layer()),
    plot_vertex_sample_(&vertex_sample_,
        history.plot().create_point(history.layer(), history.point_style()({.color_index = 15}), {vertex_sample_.w(), vertex_sample_.Lw()}),
        history.plot().create_text(history.layer(), history.label_style()({.position = cairowindow::draw::centered_below}),
          cairowindow::Point{vertex_sample_.w(), vertex_sample_.Lw()}, std::to_string(history.total_number_of_samples() - 1))) { }
};

