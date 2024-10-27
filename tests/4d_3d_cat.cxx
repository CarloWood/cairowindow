#include "sys.h"
#include "math/Polynomial.h"
#include "math/CubicPolynomial.h"
#include "gradient_descent/Sample.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Vector.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include <Eigen/Dense>
#include <thread>
#include <cmath>
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
    Window window("4d_3d_cat", 1600, 1200);

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));

    // Create another layer.
    auto second_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));

    // Open the window and start drawing.
    std::thread event_loop([&](){
      Debug(NAMESPACE_DEBUG::init_thread("event_loop"));
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    });

    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "Cubic extention of a Quartic.", {},
        "x", {},
        "y", {});
    double xz = 0;
    double yz = 0;
    double zw = 20;
    plot.set_xrange({xz - zw, xz + zw});
    plot.set_yrange({yz - zw, yz + zw});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    int color_index2 = 3; //color_pool.get_and_use_color();
    Dout(dc::notice, "color_index2 = " << color_index2);
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});
    draw::LineStyle line_style({.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::ConnectorStyle connector_style(line_style({.line_color = color::coral, .dashes = {3.0, 3.0}}));

    // Initial values.
    Point const P0{-14.8121, 13.6755}; //{-14.6654, 14.5188};
    Point const P1{-10, 15};
    Point const P2{10.5591, -14.9588};
    Point const P3{2 * P1.x() - P0.x(), 0};
    Point const P4{P2.x(), 0};

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, point_style, P0);
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, P1);
    // Create a point P₂.
    auto plot_P2 = plot.create_point(second_layer, point_style, P2);
    // Create a point P₃.
    auto plot_P3 = plot.create_point(second_layer, point_style({.filled_shape = 10}), P3);
    // Create a point P₄.
    auto plot_P4 = plot.create_point(second_layer, point_style({.filled_shape = 9}), P4);

    // Create a point Q4 on the circle.
    double const circle_radius = 1.0 * std::abs(P0.x() - P1.x());
    double phi4 = 0;
    auto plot_Q4 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P4 + circle_radius * Direction{phi4});

    // Make P0 draggable.
    window.register_draggable(plot, &plot_P0, [&plot, &plot_P1 /*, &plot_Q0*/](Point const& new_position) -> Point
        {
          Point result = new_position;
          if (new_position.x() >= plot_P1.x())
            result = Point{plot_P1.x() - 1e-6, new_position.y()};
          Dout(dc::notice, "P0 --> " << result);
          return result;
        }
    );
    // Make P1 draggable.
    window.register_draggable(plot, &plot_P1, [&plot, &plot_P0, &plot_P3 /*, &plot_Q1*/](Point const& new_position) -> Point
        {
          Point result = new_position;
          if (new_position.x() <= plot_P0.x())
            result = Point{plot_P0.x() + 1e-6, new_position.y()};
          if (new_position.x() >= plot_P3.x())
            result = Point{plot_P3.x() - 1e-6, new_position.y()};
          Dout(dc::notice, "P1 --> " << result);
          return result;
        }
    );
    // Make P2 draggable.
    window.register_draggable(plot, &plot_P2, [&](Point const& new_position) -> Point
        {
          Point result = new_position;
          if (new_position.x() <= plot_P3.x())
            result = Point{plot_P3.x() + 1e-6, new_position.y()};
          Direction const D4 = (plot_Q4 - plot_P4).direction();
          plot_P4.move(plot, {result.x(), plot_P4.y()});
          plot_Q4.move(plot, plot_P4 + circle_radius * D4);
          return result;
        }
    );
    // Make P3 draggable.
    window.register_draggable(plot, &plot_P3, [&plot, &plot_P1, &plot_P2, &plot_P3](Point const& new_position) -> Point
        {
          Point result{new_position.x(), plot_P3.y()};
          if (new_position.x() <= plot_P1.x())
            result = Point{plot_P1.x() + 1e-6, plot_P3.y()};
          if (new_position.x() >= plot_P2.x())
            result = Point{plot_P2.x() - 1e-6, plot_P3.y()};
          return result;
        }
    );
    // Make P4 draggable.
    window.register_draggable(plot, &plot_P4, [&](Point const& new_position) -> Point
        {
          //Point new_P4{plot_P2.x(), new_position.y()};
          Point new_P4{new_position};
          auto translation = new_P4 - plot_P4;
          plot_Q4.move(plot, plot_Q4 + translation);
          return new_P4;
        }
    );

    // Make Q4 draggable.
    window.register_draggable(plot, &plot_Q4, [&plot_P4, circle_radius](Point const& new_position) -> Point
        {
          Direction D4{plot_P4, new_position};
          if (D4.x() < 0.0)
            D4 = D4.inverse();
          if (D4.x() < 1e-6)
            D4 = Direction{Point{1e-6, 1.0}};
          return plot_P4 + circle_radius * D4;
        }
    );

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a label for P₀.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");

      // Draw a label for P₁.
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");

      // Draw a label for P₂.
      auto P2_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P2, "P₂");

      // Draw a circle around P₄.
      auto plot_circle4 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P4, circle_radius);

      plot::Line plot_vertical_line2{plot_P2, Direction::up};
      plot.add_line(second_layer, line_style, plot_vertical_line2);

      plot::Line plot_vertical_line3{plot_P3, Direction::up};
      plot.add_line(second_layer, line_style, plot_vertical_line3);

      // Draw a arrow from P₄ to Q₄.
      auto plot_P4_Q4 = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P4, plot_Q4);

      // Get the direction vector P₄Q₄.
      Direction const D4 = (plot_Q4 - plot_P4).direction();
      double d4 = D4.y() / D4.x();

      double x0 = plot_P0.x();
      double y0 = plot_P0.y();
      double d0 = 0.0;

      double x1 = plot_P1.x();
      double y1 = plot_P1.y();

      double x2 = plot_P2.x();
      double y2 = plot_P2.y();
      double d2 = 0.0;

      // Let f(x) be a quartic.
      math::Polynomial quartic(5 COMMA_CWDEBUG_ONLY("quartic"));

      Eigen::MatrixXd A(5, 5);

      A << std::pow(x0, 4), std::pow(x0, 3), std::pow(x0, 2), x0, 1,
         std::pow(x1, 4), std::pow(x1, 3), std::pow(x1, 2), x1, 1,
         std::pow(x2, 4), std::pow(x2, 3), std::pow(x2, 2), x2, 1,
         4 * std::pow(x0, 3), 3 * std::pow(x0, 2), 2 * x0, 1, 0,
         4 * std::pow(x2, 3), 3 * std::pow(x2, 2), 2 * x2, 1, 0;

      Eigen::VectorXd b(5);
      b << y0, y1, y2, d0, d2;

      // Step 3: Solve for the coefficients [a, b, c, d, e] using A^-1 * b
      Eigen::VectorXd coeffs = A.colPivHouseholderQr().solve(b);

      // Step 4: Output the coefficients
      for (int i = 0; i < 5; ++i)
        quartic[4 - i] = coeffs(i);

      // Plot the quartic.
//      BezierFitter fitter2([&](double w) -> Point{ return {w, quartic(w)}; }, plot.viewport());
//      auto plot_quartic = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(fitter2));

      Dout(dc::notice, "s0 = " << (Sample{plot_P0.x(), plot_P0.y(), quartic.derivative()(plot_P0.x())}) <<
          ", s1 = " << (Sample{plot_P1.x(), plot_P1.y(), quartic.derivative()(plot_P1.x())}) <<
          ", s2 = " << (Sample{plot_P2.x(), plot_P2.y(), quartic.derivative()(plot_P2.x())}));

      double a = quartic[4];
      double s = plot_P3.x() - plot_P1.x();
      math::CubicPolynomial dP(quartic.derivative());
      double d1 = dP(x1);
      math::CubicPolynomial cubic(y1, d1,
          (dP.derivative(x1) - (15.0 / 14.0) * a * s * s) / 2.0,
          (dP.second_derivative(x1) + 9.0 * a * s) / 6.0);

      // Plot the cubic.
      BezierFitter fitter3_old;
      fitter3_old.solve(
          [&, x1](double w) -> Point { return {w, cubic(w - x1)}; },
          plot.viewport());
      auto old_plot_cubic = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::green}), std::move(fitter3_old));

      BezierFitter fitter3;
      Rectangle viewport{-20, -20, 40, 40};
      fitter3.solve(
          [&, x1](double w) -> Point { return {w, cubic(w - x1)}; },
          [&, x1](double w) -> Vector { return {1.0, cubic.derivative(w - x1)}; },
          {viewport.offset_x(), viewport.offset_x() + viewport.width()},
          viewport, 1e-5, Orientation::vertical);
      auto plot_cubic = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::red}), std::move(fitter3));

      math::CubicPolynomial cubic2;
      cubic2.initialize(plot_P1.x(), plot_P1.y(), d1, plot_P4.x(), plot_P4.y(), d4);

      // Plot the cubic2.
//      BezierFitter fitter4([&](double w) -> Point{ return {w, cubic2(w)}; }, plot.viewport());
//      auto plot_cubic2 = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::green}), std::move(fitter4));

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
