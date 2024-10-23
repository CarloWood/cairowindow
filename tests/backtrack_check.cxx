#include "sys.h"
#include "math/Polynomial.h"
#include "math/QuadraticPolynomial.h"
#include "cairowindow/BezierCurve.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
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
    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "Backtrack Check", {},
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

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, point_style, {-6.0, 4.0});
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, {2.0, 0.0});

    // Create a point Q on the circle.
    double phi = -0.25 * M_PI;
    double const circle_radius = 2.0;
    auto plot_Q = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P0 + circle_radius * Direction{phi});

    // Make all points draggable.
    window.register_draggable(plot, &plot_P0);
    window.register_draggable(plot, &plot_P1, [&plot_P1](Point const& new_position) -> Point
        {
          return {plot_P1.x(), new_position.y()};
        }
    );

    window.register_draggable(plot, &plot_Q, [&plot_P0, circle_radius](Point const& new_position)
        {
          Direction D0{plot_P0, new_position};
          if (D0.y() > 0.0)
            D0 = D0.inverse();
          if (D0.y() > -10e-6)
            D0 = Direction{Point{1.0, -10e-6}};
          return plot_P0 + circle_radius * D0;
        }
    );

    auto slider_c = plot.create_slider(second_layer, {978, 83, 7, 400}, 0.3, 0.01, 1.5);
    auto slider_c_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "c");

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a circle around P₀.
      auto plot_circle = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P0, circle_radius);

      // Draw a label for P₀.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");

      // Draw a arrow from P₀ to Q.
      auto plot_P0_Q = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P0, plot_Q);

      // Store the velocity vectors as Vector.
      Direction const D0 = (plot_Q - plot_P0).direction();

      // Get w₀.
      double w0 = plot_P0.x();
      // Get L(w₀)
      double Lw0 = plot_P0.y();
      // Get L'(w₀)
      double dLdw0 = D0.y() / D0.x();
      // Get c.
      double c = slider_c.value();

      // P(w)    = a + b w + c w²
      // P(w₀)   = a + b w₀ + c w₀² = L(w₀)     -->  a = L(w₀) - b w₀ - c w₀²
      // P'(w₀)  =     b    + 2c w₀ = L'(w₀)    -->  b = L'(w₀) - 2c w₀
      double b = dLdw0 - 2.0 * c * w0;
      double a = Lw0 - b * w0 - c * utils::square(w0);
      math::QuadraticPolynomial parabola{a, b, c};

      // Plot the parabola.
      BezierFitter fitter([&](double w) -> Point{ return {w, parabola(w)}; }, plot.viewport());
      auto plot_parabola = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(fitter));

      // Draw a vertical line through the vertex of the parabola.
      auto plot_vertex_line = plot.create_line(second_layer, line_style, Point{parabola.vertex_x(), 0.0}, Direction::up);

      // Force P1 on that line.
      plot_P1.move(plot, Point{parabola.vertex_x(), plot_P1.y()});
      Dout(dc::notice, "parabola.vertex_x() = " << parabola.vertex_x() << ", plot_P1 = " << plot_P1);

      // Draw a label for P₁.
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");

      // Draw line through P₀ and P₁.
      auto plot_line_P0P1 = plot.create_line(second_layer, line_style, plot::LineExtend::both, plot_P0, plot_P1);

      // Get w₁.
      double w1 = plot_P1.x();
      // Get L(w₁)
      double Lw1 = plot_P1.y();

      Dout(dc::notice, "w₀ = " << w0 << "; L(w₀) = " << Lw0 << "; L'(w₀) = " << dLdw0 <<
          "; w₁ = " << w1 << "; L(w₁) = " << Lw1 << "; c = " << c);

      // Let Q(w) = d + e w + f w² + g w³ be a cubic that has the same value and first and
      // second derivative in w₀, as well as goes through the point (w₁, L(w₁)).
      math::Polynomial cubic(4 COMMA_CWDEBUG_ONLY("cubic"));
      // Then (see README.gradient_descent),
      //
      //      L(w₀) - L(w₁) - (w₀ - w₁) L'(w₀) + (w₀ - w₁)² c
      // g =  -----------------------------------------------
      //                      (w₀ - w₁)³
      //
      // f = c - 3g w₀
      //
      // e = L'(w₀) - 2f w₀ - 3g w₀²
      //
      double dw = w0 - w1;
      double dw2 = utils::square(dw);
      double dw3 = dw * dw2;
      double g = (Lw0 - Lw1 - dw * dLdw0 + dw2 * c) / dw3;
      double f = c - 3.0 * g * w0;
      double e = dLdw0 - 2.0 * f * w0 - 3.0 * g * utils::square(w0);

      Dout(dc::notice, "g = " << g);

      cubic[3] = g;
      cubic[2] = f;
      cubic[1] = e;
      cubic[0] = Lw0 - cubic(w0);

      double h = Lw1 - parabola(w1);
      double alpha = std::pow(2.0 * c / (b + 2 * c * w0), 3.0);

      // Plot the cubic.
      BezierFitter fitter2([&](double w) -> Point{ return {w, /*cubic(w)*/ parabola(w) - h * alpha * std::pow(w - w0, 3.0)}; }, plot.viewport());
      auto plot_cubic = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(fitter2));

      double beta = 6.0 * h * alpha;
//      double zero = (2.0 * c + beta * w0 - std::sqrt(utils::square(2.0 * c) + 2.0 * beta * (b + 2.0 * c * w0))) / beta;

      double gamma = beta * (w0 - w1) / c;
      double zero = w0 + 2.0 * (std::sqrt(1.0 + gamma) - 1) / gamma * (w1 - w0);

      // Draw a vertical line where the zero is.
      auto plot_zero_line = plot.create_line(second_layer, line_style, Point{zero, 0.0}, Direction::up);

      Dout(dc::notice, "gamma = " << gamma);

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
