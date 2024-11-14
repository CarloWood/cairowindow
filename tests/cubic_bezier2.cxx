#include "sys.h"
#include "math/Polynomial.h"
#include "cairowindow/BezierCurve.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Matrix.h"
#include "cairowindow/draw/Shape.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/draw/Point.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include "utils/square.h"
#include "utils/almost_equal.h"
#include <Eigen/Dense>
#include <thread>
#include <iostream>
#include "debug.h"

int main()
{
#ifdef CWDEBUG
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Cubic Bezier", 1200, 900);

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
        "Cubic Bezier test", {},
        "x", {},
        "y", {});
    plot.set_xrange({-10, 10});
    plot.set_yrange({-10, 10});
    plot.add_to(background_layer, true);

    // Create a second window.
    Window window_newton_raphson("Newton-Raphson", 1200, 900);

    // Create a new layer with a gray background.
    auto background_layer_newton_raphson =
      window_newton_raphson.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer_newton_raphson"));

    // Create another layer.
    auto second_layer_newton_raphson = window_newton_raphson.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer_newton_raphson"));

    // Open the window and start drawing.
    std::thread event_loop_newton_raphson([&](){
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window_newton_raphson.run();
      event_loop.set_cleanly_terminated();
    });

    // Create and draw plot area.
    plot::Plot plot_newton_raphson(window_newton_raphson.geometry(), { .grid = {.color = color::orange} },
        "Newton Raphson test", {},
        "g", {},
        "dJ/dg", {});
    plot_newton_raphson.set_xrange({0.0, 1.0});
    plot_newton_raphson.set_yrange({-1e6, 1e6});
    plot_newton_raphson.add_to(background_layer_newton_raphson, false);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    draw::PointStyle C0_point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 12});
    draw::PointStyle C1_point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 12});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});
    draw::LineStyle line_style({.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::ArcStyle arc_style({.line_color = color::blue, .line_width = 1.0});
    draw::RectangleStyle rectangle_style({.line_color = color::red, .line_width = 1.0});

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, point_style, {-9.50695, 0.0505689});
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, {2.88243, -3.43869});
    // Create a point Pᵧ.
    auto plot_P_gamma = plot.create_point(second_layer, point_style, {-8.87484, 1.03666});
    // Create a point Q.
    auto plot_Q = plot.create_point(second_layer, point_style, {0.151707, 1.64349});

    Vector initial_T0{1.07732, 1.68504};
    Vector initial_T1{0.218681, 1.98801};

    // Create point Q0 and Q1 on circles around P0 and P1.
    double phi0 = initial_T0.direction().as_angle();
    double phi1 = initial_T1.direction().as_angle();

    double const circle_radius = 2.0;
    auto plot_Q0 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P0 + circle_radius * Direction{phi0});
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
    window.register_draggable(plot, &plot_P_gamma);
    window.register_draggable(plot, &plot_Q);

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

//    auto slider_gamma = plot.create_slider(second_layer, {978, 83, 7, 400}, 0.5, 0.01, 0.99);
//    auto slider_gamma_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "γ");

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);
      window_newton_raphson.set_send_expose_events(false);

      // Draw a circle around P₀.
      auto plot_circle0 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P0, circle_radius);

      // Draw a label for P₀.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");

      // Draw a arrow from P₀ to Q.
      auto plot_P0_Q = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P0, plot_Q0);

      // Draw a circle around P₁.
      auto plot_circle1 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P1, circle_radius);

      // Draw a label for P₁.
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");

      // Draw a arrow from P₁ to Q.
      auto plot_P1_Q = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P1, plot_Q1);

      // Draw a label for Pᵧ.
      auto P_gamma_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P_gamma, "Pᵧ");

      // Draw a label for Q.
      auto P_Q_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_Q, "Q");

      Dout(dc::notice, "P0 = " << plot_P0 << "; P1 = " << plot_P1 << "; P_gamma = " << plot_P_gamma);

      // Get the distance between P0 and P1.
      double dist = (plot_P1 - plot_P0).length();

      // Calculate the required tangent vectors, normalized wrt the distance between P0 and P1.
      Vector const T0 = plot_Q0 - plot_P0;
      Vector const T1 = plot_Q1 - plot_P1;

      Dout(dc::notice, "T0 = " << T0 << "; T1 = " << T1);

      BezierCurve bezier_curve(plot_P0, plot_P1);
      Vector S = bezier_curve.cubic_from(T0, T1, plot_P_gamma);
      double const& Sx = S.x();
      double const& Sy = S.y();

      Dout(dc::notice, "α T0 = " << (Sx * T0) << "; β T1 = " << (Sy * T1));

      bool reject = Sx < 0.1 || Sy < 0.1;

      // Draw the Bezier curve.
      //plot_bezier_curve = plot.create_bezier_curve(second_layer, curve_line_style, bezier_curve);

      plot::BezierFitter plot_bezier_fitter;

      // The life-time this object determines how long the curve is visible.
      BezierFitter bezier_fitter(
          [&bezier_curve](double t) -> Point
          {
            //Point p{3.0 * std::cos(2.0 * M_PI * t), 3.0 * std::sin(2.0 * M_PI * t)};
            //Dout(dc::notice, "P(" << t << ") = " << p);
            //return p;
            return bezier_curve.P(t);
          },
          {0.0, 1.0},           // Domain
          plot.viewport());

      Dout(dc::notice, "bezier_curve = " << bezier_curve);
      Dout(dc::notice, "|J| = " << bezier_curve.J().length());

      // Draw bezier_fitter.
      plot_bezier_fitter = plot.create_bezier_fitter(second_layer,
          curve_line_style({.line_color = reject ? color::red : color::black}), std::move(bezier_fitter));

      plot::Point plot_C0 = plot.create_point(second_layer, C0_point_style, bezier_curve.C0().as_point());
      plot::Point plot_C1 = plot.create_point(second_layer, C1_point_style, bezier_curve.C1().as_point());

      BezierCurve quadratic_bezier_curve(plot_P0, plot_P1);
      quadratic_bezier_curve.quadratic_from(T0.direction(), T1.direction());
      auto plot_quadratic_bezier_curve = plot.create_bezier_curve(second_layer, line_style, quadratic_bezier_curve);

      plot::BezierFitter plot_dJdg;
      std::unique_ptr<CubicBezierCurve> cubic_bezier_curve = CubicBezierCurve::create(plot_P0, plot_P1, T0, T1, plot_P_gamma);
      BezierFitter bezier_fitter_dJdg(
          [&](double g) -> Point
          {
            return {g, cubic_bezier_curve->derivative(g)};
          },
          {0.01, 0.99},
          plot_newton_raphson.viewport());
      plot_dJdg = plot_newton_raphson.create_bezier_fitter(second_layer_newton_raphson,
          curve_line_style, std::move(bezier_fitter_dJdg));

      plot::BezierFitter plot_ddJdgdg;
      BezierFitter bezier_fitter_ddJdgdg(
          [&](double g) -> Point
          {
            return {g, 0.1 * cubic_bezier_curve->second_derivative(g)};
          },
          {0.01, 0.99},
          plot_newton_raphson.viewport());
      plot_ddJdgdg = plot_newton_raphson.create_bezier_fitter(second_layer_newton_raphson,
          curve_line_style({.line_color = color::red}), std::move(bezier_fitter_ddJdgdg));

      cubic_bezier_curve.reset();

#if 0
      Debug(libcw_do.off());
      // Draw J.
      BezierFitter J(
          [&](double g) -> Point
          {
            BezierCurve bc(plot_P0, plot_P1);
            [[maybe_unused]] Vector S = bc.cubic_from(T0, T1, plot_P_gamma, g);
            return (0.01 * bc.J()).as_point();
          },
          {0.01, 0.99},
          plot.viewport()
          );
      plot::BezierFitter plot_J_curve = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::blue}), std::move(J));
      Debug(libcw_do.on());
#endif

#if 1
      Vector const J = bezier_curve.J();
      Vector const A0 = bezier_curve.A0();
      Vector const V0 = bezier_curve.V0();
      Vector const B{plot_P0};        // = bezier_curve.B();      // This is the Beginning of the curve, aka P0.

      // Solve P'⋅(P - Q) = 0, to find all t such that P(t) is a point on the Bezier curve
      // where the tangent of the curve is perpendicular to a line from Q to that point P.
      //
      //  P(t) = B + V0 t + (A0/2) t² + J/6 t³
      // P'(t) = V0 + A0 t + J/2 t²
      //
      // P - Q = (B - Q) + V0 t + (A0/2) t² + J/6 t³

      // Abbreviate B - Q.
      Vector const QB = plot_P0 - plot_Q;

      //
      // 0 = P'⋅(P - Q) =     V0⋅QB    +    V0⋅V0 t  + 1/2 V0⋅A0 t² + 1/6 V0⋅J t³ +
      //                      A0⋅QB t  +    A0⋅V0 t² + 1/2 A0⋅A0 t³ + 1/6 A0⋅J t⁴ +
      //                   1/2 J⋅QB t² + 1/2 J⋅V0 t³ + 1/4  J⋅A0 t⁴ + 1/12 J⋅J t⁵
      //
      // Multiply with 12 to get:
      //
      // 0 =               12 V0⋅QB    + 12 V0⋅V0 t  +   6 V0⋅A0 t² +   2 V0⋅J t³ +
      //                   12 A0⋅QB t  + 12 A0⋅V0 t² +   6 A0⋅A0 t³ +   2 A0⋅J t⁴ +
      //                    6  J⋅QB t² +  6  J⋅V0 t³ +   3  J⋅A0 t⁴ +      J⋅J t⁵
      //
      // Thus,
      //
      // 0 = J² t⁵ + 5 J⋅A0 t⁴ + (8 J⋅V0 + 6 A0²) t³ + (6 J⋅QB + 18 A0⋅V0) t² + (12 A0⋅QB + 12 V0²) t + 12 V0⋅QB

      math::Polynomial polynomial(6 COMMA_CWDEBUG_ONLY("t"));
      polynomial[5] = J.length_squared();
      polynomial[4] = 5.0 * J.dot(A0);
      polynomial[3] = 8.0 * J.dot(V0) + 6.0 * A0.length_squared();
      polynomial[2] = 6.0 * (J.dot(QB) + 3.0 * A0.dot(V0));
      polynomial[1] = 12.0 * (A0.dot(QB) + V0.length_squared());
      polynomial[0] = 12.0 * V0.dot(QB);

      // Get the real roots of polynomial to find the points Pt on the curve where PtQ makes an angle of 90 degrees with the curve.
      std::array<std::complex<double>, 5> roots;
      int number_of_roots = polynomial.get_roots(roots);

      // Find the shorest distance.
      double min_distance_squared = std::numeric_limits<double>::max();
      int min_r = -1;
      for (int r = 0; r < number_of_roots; ++r)
      {
        bool is_real = utils::almost_equal(roots[r], std::conj(roots[r]), 1e-4);
        if (!is_real)
          continue;
        double t = roots[r].real();
        if (t < 0.0 || t > 1.0)
          continue;
        double distance_squared = (bezier_curve.P(t) - plot_Q).length_squared();
        if (distance_squared < min_distance_squared)
        {
          min_distance_squared = distance_squared;
          min_r = r;
        }
      }
      // Also calculate the distance to P₀ and P₁.
      double distance_P0_squared = QB.length_squared();
      double distance_P1_squared = (plot_P1 - plot_Q).length_squared();
      // If one of those distances is less than the shortest distance found so far, than
      // this point Q does not belong to this Bezier curve.
      plot::Connector plot_shortest_distance;
      if (distance_P0_squared > min_distance_squared && distance_P1_squared > min_distance_squared)
      {
        // Draw a line for the shortest distance found, from Q to P(root[min_r]).
        plot_shortest_distance = plot.create_connector(second_layer, solid_line_style({.line_color = color::lime}),
            plot_Q, bezier_curve.P(roots[min_r].real()));
      }
#endif

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);
      window_newton_raphson.set_send_expose_events(true);

      // Block until a redraw is necessary (for example because the user moved a draggable object,
      // or wants to print the current drawing) then go to the top of loop for a redraw.
      if (!window.handle_input_events())
      {
        window_newton_raphson.close();
        break;          // Program must be terminated.
      }
    }

    event_loop.join();
    event_loop_newton_raphson.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
#else
  std::cerr << "Requires Debug mode." << std::endl;
#endif
}
