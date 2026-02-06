#include "sys.h"
#include "math/QuadraticPolynomial.h"
#include "math/Polynomial.h"
#include "cairowindow/BezierCurve.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Direction.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/draw/Point.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include "utils/square.h"
#include <thread>
#include "debug.h"

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Backtrack Check", 1200, 900);

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));

    // Create another layer.
    auto second_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));

    // Open the window and start drawing.
    std::thread event_loop([&](){
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    });

    // Create and draw plot area.
    plot::Plot plot(window.geometry(), draw::PlotAreaStyle({.color = color::orange}),
        "Parabolic fit", {},
        "x", {},
        "y", {});
    plot.set_xrange({-10, 10});
    plot.set_yrange({-10, 10});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});
    draw::LineStyle line_style({.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::ArcStyle arc_style({.line_color = color::blue, .line_width = 1.0});
    draw::ConnectorStyle error_style({.line_color = color::red});

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, point_style, {-6.0, 4.0});
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, {3.0847, -3.51454});

    // Create a point Q0 on the circle.
//    double phi0 = -0.25 * M_PI;
    double phi0 = -1.19007;

    double const circle_radius = 2.0;
    auto plot_Q0 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P0 + circle_radius * Direction{phi0});
    // Create a point Q1 on the circle.
//    double phi1 = 0.1 * M_PI;
    double phi1 = -0.234898;
    auto plot_Q1 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P1 + circle_radius * Direction{phi1});

    // Make all points draggable.
    window.register_draggable(plot, &plot_P0, [&plot, &plot_P0, &plot_Q0](Point const& new_position) -> Point
        {
          auto translation = new_position - plot_P0;
          plot_Q0.move(plot, plot_Q0 + translation);
          return new_position;
        }
    );
    window.register_draggable(plot, &plot_P1, [&plot, &plot_P1, &plot_Q1](Point const& new_position) -> Point
        {
          auto translation = new_position - plot_P1;
          plot_Q1.move(plot, plot_Q1 + translation);
          return new_position;
        }
    );

    window.register_draggable(plot, &plot_Q0, [&plot_P0, circle_radius](Point const& new_position) -> Point
        {
          return plot_P0 + circle_radius * Direction{plot_P0, new_position};
        }
    );
    window.register_draggable(plot, &plot_Q1, [&plot_P1, circle_radius](Point const& new_position) -> Point
        {
          return plot_P1 + circle_radius * Direction{plot_P1, new_position};
        }
    );

    auto slider_c = plot.create_slider(second_layer, {978, 83, 7, 400}, 0.0520743, -1.0, 1.0);
    auto slider_c_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "c");

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a circle around P₀.
      auto plot_circle0 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P0, circle_radius);

      // Draw a label for P₀.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");

      // Draw a arrow from P₀ to Q.
      auto plot_P0_Q = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P0, plot_Q0);

      // Store the velocity vectors as Vector.
      Direction const D0 = (plot_Q0 - plot_P0).direction();

      Dout(dc::notice, "phi0 = " << D0.as_angle());

      // Draw a circle around P₁.
      auto plot_circle1 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P1, circle_radius);

      // Draw a label for P₁.
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");

      // Draw a arrow from P₁ to Q.
      auto plot_P1_Q = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P1, plot_Q1);

      // Store the velocity vectors as Vector.
      Direction const D1 = (plot_Q1 - plot_P1).direction();

      Dout(dc::notice, "phi1 = " << D1.as_angle());

      // Draw line through P₀ and P₁.
      auto plot_line_P0P1 = plot.create_line(second_layer, line_style, LineExtend::both, plot_P0, plot_P1);

      // Get w₀.
      double w0 = plot_P0.x();
      // Get L(w₀)
      double Lw0 = plot_P0.y();
      // Get L'(w₀)
      double dLdw0 = D0.y() / D0.x();
      // Get w₁.
      double w1 = plot_P1.x();
      // Get L(w₁)
      double Lw1 = plot_P1.y();
      // Get L'(w₁)
      double dLdw1 = D1.y() / D1.x();

      double delta_w = w0 - w1;
      double delta_L = Lw0 - Lw1;
      double delta_dLdw = dLdw0 - dLdw1;
      double sum_dLdw = dLdw0 + dLdw1;
      double product_L = dLdw0 * dLdw1;

      // Get c.
      double c_s = slider_c.value();
      double c_a = product_L < 0.0 ?
        (delta_L * sum_dLdw - 2.0 * delta_w * product_L) / (utils::square(delta_w) * delta_dLdw) :
        delta_L * delta_dLdw / (utils::square(delta_w) * sum_dLdw);

      Dout(dc::notice, "w₀ = " << w0 << "; L(w₀) = " << Lw0 << "; L'(w₀) = " << dLdw0 <<
          "; w₁ = " << w1 << "; L(w₁) = " << Lw1 << "; c = " << c_a);

      double b_s = (Lw0 - Lw1) / (w0 - w1) - (w0 + w1) * c_s;
      double b_a = (Lw0 - Lw1) / (w0 - w1) - (w0 + w1) * c_a;
      double a_s = (w1 * Lw0 - w0 * Lw1) / (w1 - w0) + w0 * w1 * c_s;
      double a_a = (w1 * Lw0 - w0 * Lw1) / (w1 - w0) + w0 * w1 * c_a;
      math::QuadraticPolynomial<double> parabola_s{a_s, b_s, c_s};
      math::QuadraticPolynomial<double> parabola_a{a_a, b_a, c_a};

      // Plot the parabola.
      BezierFitter fitter_s([&](double w) -> Point{ return {w, parabola_s(w)}; }, plot.viewport());
      auto plot_parabola_s = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(fitter_s));
      // Plot the parabola.
      BezierFitter fitter_a([&](double w) -> Point{ return {w, parabola_a(w)}; }, plot.viewport());
      auto plot_parabola_a = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(fitter_a));

      // Draw a vertical line through the vertex of the parabola_a.
      auto plot_vertex_line = plot.create_line(second_layer, line_style, Point{parabola_a.vertex_x(), 0.0}, Direction::up);

      // Draw lines tangent to parabola_a.
      auto plot_tangent_P0 = plot.create_line(second_layer, solid_line_style({.line_color = color::lightgray}),
          plot_P0, Direction{std::atan(parabola_a.derivative(w0))});
      auto plot_tangent_P1 = plot.create_line(second_layer, solid_line_style({.line_color = color::lightgray}),
          plot_P1, Direction{std::atan(parabola_a.derivative(w1))});

      // Calculate the relative error of the derivative at w₀.
      double re0_s = (parabola_s.derivative(w0) - dLdw0) / dLdw0;
      // Calculate the relative error of the derivative at w₁.
      double re1_s = (parabola_s.derivative(w1) - dLdw1) / dLdw1;
      // Calculate the relative error of the derivative at w₀.
      double re0_a = (parabola_a.derivative(w0) - dLdw0) / dLdw0;
      // Calculate the relative error of the derivative at w₁.
      double re1_a = (parabola_a.derivative(w1) - dLdw1) / dLdw1;

      // Draw an arrow representing the relative errors.
      plot::Connector error0_a(plot_P0, plot_P0 + Vector{0.0, 5.0 * re0_a});
      if (!std::isnan(re0_a) && !std::isinf(re0_a))
        plot.add_connector(second_layer, error_style, error0_a);
      plot::Connector error1_a(plot_P1, plot_P1 + Vector{0.0, 5.0 * re1_a});
      if (!std::isnan(re1_a) && !std::isinf(re1_a))
        plot.add_connector(second_layer, error_style, error1_a);
      // Draw an arrow representing the relative errors.
      plot::Connector error0_s(plot_P0, plot_P0 + Vector{0.0, 5.0 * re0_s});
      if (!std::isnan(re0_s) && !std::isinf(re0_s))
        plot.add_connector(second_layer, error_style({.line_color = color::blue, .dashes = {10.0, 10.0}}), error0_s);
      plot::Connector error1_s(plot_P1, plot_P1 + Vector{0.0, 5.0 * re1_s});
      if (!std::isnan(re1_s) && !std::isinf(re1_s))
        plot.add_connector(second_layer, error_style({.line_color = color::blue, .dashes = {10.0, 10.0}}), error1_s);

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // Block until a redraw is necessary (for example because the user moved a draggable object,
      // or wants to print the current drawing) then go to the top of loop for a redraw.
      if (!window.handle_input_events())
        break;          // Program must be terminated.
    }

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
