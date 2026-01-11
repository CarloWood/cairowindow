#include "sys.h"
#include "math/QuadraticPolynomial.h"
#include "math/Polynomial.h"
#include "gradient_descent/Sample.h"
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
    using Sample = gradient_descent::Sample;

    // Create a window.
    Window window("Clamp Check", 1200, 900);

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

    Sample s1{10.1586, 47.5865, -0.164047}, s2{11.8608, 47.013, -0.358701};

    double w1 = std::min(s1.w(), s2.w());
    double w2 = std::max(s1.w(), s2.w());
    double dw = w2 - w1;
    double const w_min = w1 - 0.2 * dw;
    double const w_max = w2 + 0.2 * dw;

    double Lw1 = std::min(s1.Lw(), s2.Lw());
    double Lw2 = std::max(s1.Lw(), s2.Lw());
    double dLw = Lw2 - Lw1;
    double const Lw_min = Lw1 - 0.2 * dLw;
    double const Lw_max = Lw2 + 0.2 * dLw;

    // Create and draw plot area.
    plot::Plot plot(window.geometry(), draw::PlotAreaStyle({.color = color::orange}),
        "Clamp Check", {},
        "x", {},
        "y", {});
    plot.set_xrange({w_min, w_max});
    plot.set_yrange({Lw_min, Lw_max});
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
    auto plot_P0 = plot.create_point(second_layer, point_style, {s1.w(), s1.Lw()});
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, {s2.w(), s2.Lw()});

    // Create a point Q0 on the circle.
    double const circle_radius = 0.2 * std::abs(s2.w() - s1.w());
    double phi0 = std::atan(s1.dLdw());
    auto plot_Q0 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P0 + circle_radius * Direction{phi0});
    double phi1 = std::atan(s2.dLdw());
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
          Direction D0{plot_P0, new_position};
          if (D0.x() < 0.0)
            D0.negate();
          if (D0.x() < 1e-6)
            D0 = Direction{Point{1e-6, 1.0}};
          return plot_P0 + circle_radius * D0;
        }
    );
    window.register_draggable(plot, &plot_Q1, [&plot_P1, circle_radius](Point const& new_position) -> Point
        {
          Direction D1{plot_P1, new_position};
          if (D1.x() < 0.0)
            D1.negate();
          if (D1.x() < 1e-6)
            D1 = Direction{Point{1e-6, 1.0}};
          return plot_P1 + circle_radius * D1;
        }
    );

//    auto slider_d = plot.create_slider(second_layer, {978, 83, 7, 400}, 0.0520743, -0.01, 0.01);
//    auto slider_d_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "d");

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a circle around P₀.
      auto plot_circle0 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P0, circle_radius);

      // Draw a label for P₀.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");

      // Draw a arrow from P₀ to Q₀.
      auto plot_P0_Q0 = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P0, plot_Q0);

      // Get the direction vector of L'(w₀).
      Direction const D0 = (plot_Q0 - plot_P0).direction();

      // Draw a circle around P₁.
      auto plot_circle1 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P1, circle_radius);

      // Draw a label for P₁.
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");

      // Draw a arrow from P₁ to Q₁.
      auto plot_P1_Q1 = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P1, plot_Q1);

      // Get the direction vector of L'(w₁).
      Direction const D1 = (plot_Q1 - plot_P1).direction();

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

      Dout(dc::notice, "w₀ = " << w0 << "; L(w₀) = " << Lw0 << "; L'(w₀) = " << dLdw0 <<
          "; w₁ = " << w1 << "; L(w₁) = " << Lw1 << "; L'(w₁) = " << dLdw1);

      // Let Q(w) = a + b w + c w² + d w³ be a cubic that has the same value and first derivative in w₀ and w₁.
      math::Polynomial<double> cubic(4 COMMA_CWDEBUG_ONLY("cubic"));

      // See https://math.stackexchange.com/questions/4926335/
      double dw = w0 - w1;
      double dw3 = std::pow(dw, 3.0);
      double d = (-2.0 * (Lw0 - Lw1) + dw * (dLdw0 + dLdw1)) / dw3;
      double c = (dLdw0 - dLdw1) / (2.0 * dw) - 1.5 * (w0 + w1) * d;
      double b = (w0 * dLdw1 - w1 * dLdw0) / dw + 3.0 * w0 * w1 * d;

      cubic[3] = d;
      cubic[2] = c;
      cubic[1] = b;
      cubic[0] = Lw0 - cubic(w0);

      // Plot the cubic.
      BezierFitter fitter2([&](double w) -> Point{ return {w, cubic(w)}; }, plot.viewport());
      auto plot_cubic = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(fitter2));

      double D = utils::square(c) - 3.0 * b * d;

      plot::Line plot_minimum_line;
      plot::Line plot_maximum_line;
      if (D >= 0.0)
      {
        // Find the minimum (the maximum would be `- std::sqrt`).
        double minimum = (-c + std::sqrt(D)) / (3.0 * d);
        // Draw a vertical line where the minimum is.
        plot_minimum_line = plot.create_line(second_layer, line_style, Point{minimum, 0.0}, Direction::up);

        double maximum = (-c - std::sqrt(D)) / (3.0 * d);
        // Draw a vertical line where the maximum is.
        plot_maximum_line = plot.create_line(second_layer, line_style, Point{maximum, 0.0}, Direction::up);

        Dout(dc::notice, std::setprecision(std::numeric_limits<long double>::max_digits10) <<
            "minimum = " << minimum << ", maximum = " << maximum);
      }

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
