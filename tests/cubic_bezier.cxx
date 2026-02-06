#include "sys.h"
#include "cairowindow/BezierCurve.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Direction.h"
#include "cairowindow/draw/Shape.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/draw/Point.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include "utils/square.h"
#include <thread>
#include <iostream>
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
    plot::Plot plot(window.geometry(), draw::PlotAreaStyle({.color = color::orange}),
        "Cubic Bezier test", {},
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
    draw::RectangleStyle rectangle_style({.line_color = color::red, .line_width = 1.0});

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, point_style, {-6.0, -1.0});
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, {2.0, 0.0});

    // Create a point D₀.
//    auto plot_D0 = plot.create_point(second_layer, point_style({.color_index = 2}), {2.5, 1.5});
    // Create a point D₁.
//    auto plot_D1 = plot.create_point(second_layer, point_style({.color_index = 2}), {3.5, -1.5});

    // Create a point Pᵧ.
    plot::Point plot_P_gamma = plot.create_point(second_layer, point_style({.filled_shape = 10}), {-4.0, 1.0});

    // Create a point Q on the circle.
    double phi = 0.3 * M_PI;
    double const circle_radius = 2.0;
    auto plot_Q = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P_gamma + circle_radius * Direction{phi});

    // Make all points draggable.
    window.register_draggable(plot, &plot_P0);
    window.register_draggable(plot, &plot_P1);
//    window.register_draggable(plot, &plot_D0);
//    window.register_draggable(plot, &plot_D1);
    window.register_draggable(plot, &plot_P_gamma);

    window.register_draggable(plot, &plot_Q, [&plot_P_gamma, circle_radius](Point const& new_position) -> Point
        {
          return plot_P_gamma + circle_radius * Direction{plot_P_gamma, new_position};
        }
    );

//    auto slider_k0 = plot.create_slider(second_layer, {928, 83, 7, 400}, 0.25, -10.0, 10.0);
//    auto slider_k0_label = plot.create_text(second_layer, slider_style, Pixel{928, 483}, "k0");

    auto slider_beta = plot.create_slider(second_layer, {978, 83, 7, 400}, 0.0, -100.0, 100.0);
    auto slider_beta_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "beta");

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a circle around Pᵧ.
      auto plot_circle = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P_gamma, circle_radius);

      // Draw a label for P₀, P₁, D₀ and D₁.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");
//      auto D0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_D0, "D₀");
//      auto D1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_D1, "D₁");
      auto P_gamma_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P_gamma, "Pᵧ");

      // Draw a arrow from Pᵧ to Q.
      auto plot_Pgamma_Q = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P_gamma, plot_Q);
      // Draw a label for Q.
      auto Q_label = plot.create_text(second_layer, label_style, plot_Q, "Q");

      // Draw line through P₀ and Pᵧ.
      auto plot_line_P0Pgamm = plot.create_line(second_layer, line_style, LineExtend::both, plot_P0, plot_P_gamma);

      // Draw line through P₁ and Pᵧ.
      auto plot_line_P1Pgamma = plot.create_line(second_layer, line_style, LineExtend::both, plot_P1, plot_P_gamma);

      // Draw an arrow from P₀ to D₀.
//      auto plot_D0_arrow = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P0, plot_D0);
      // Draw an arrow from P₁ to D₁.
//      auto plot_D1_arrow = plot.create_connector(second_layer, line_style({.line_color = color::orange}), plot_P1, plot_D1);

      // Store the velocity vectors as Vector.
      Direction const D_gamma = (plot_Q - plot_P_gamma).direction();
//      Vector const V1 = plot_D1_arrow;

      plot::Connector plot_curvature;
      Vector K0(0.0, 0.0);
      Rectangle rect;

#if 0
      std::vector<Point> curve_points;
      {
        // Define the matrix.

#if 1
        // Using α·D₀ and β·D₁.

        // X(t) = P₀ + α·D₀·t + (3(P₁-P₀)-2α·D₀-β·D₁)·t² + (-2(P₁-P₀)+α·D₀+β·D₁)·t³
        //
        //      ⎡Px₀  Vx₀  (3(Px₁-Px₀)-2Vx₀-Vx₁)  (-2(Px₁-Px₀)+Vx₀+Vx₁)⎤
        //  M = ⎣Py₀  Vy₀  (3(Py₁-Py₀)-2Vy₀-Vy₁)  (-2(Py₁-Py₀)+Vy₀+Vy₁)⎦
        //

        Vector V2 = 3 * (plot_P1 - plot_P0) - 2 * V0 - V1;
        Vector V3 = -2 * (plot_P1 - plot_P0) + V0 + V1;
#else
        // Using α·D₀ and k₀.

        //     ⎡Px₀  α·Dx₀  (β·Dx₀ - α²k₀·Dy₀)/2  Px₁-(Px₀+α·Dx₀+(β·Dx₀ - α²k₀·Dy₀)/2)⎤
        // M = ⎣Py₀  α·Dy₀  (β·Dy₀ + α²k₀·Dx₀)/2  Py₁-(Py₀+α·Dy₀+(β·Dy₀ + α²k₀·Dx₀)/2)⎦

        double alpha = V0.norm();
        Direction D0 = V0.direction();
        Vector V2 = 0.5 * (slider_beta.value() * D0 + alpha * alpha * slider_k0.value() * D0.normal());
        Vector V3 = plot_P1 - plot_P0 - V0 - V2;
#endif

        double m00 = plot_P0.x();
        double m10 = plot_P0.y();
        double m01 = V0.x();
        double m11 = V0.y();
        double m02 = V2.x();
        double m12 = V2.y();
        double m03 = V3.x();
        double m13 = V3.y();

        BezierCurveMatrix M = {{{{m00, m10}, {m01, m11}, {m02, m12}, {m03, m13}}}};
        BezierCurve bc(M);
        rect = bc.extents();

        auto xt = [=](double t){ return m00 + t * (m01 + t * (m02 + t * m03)); };
        auto yt = [=](double t){ return m10 + t * (m11 + t * (m12 + t * m13)); };

        for (int i = -200; i <= 400; ++i)
        {
          double t = i * 0.01;
          curve_points.emplace_back(xt(t), yt(t));
        }

        // Now we have the matrix, lets calculate the second derivative of X(t).
        // Zero's derivative:
        // x(t) = m00 + m01 t + m02 t^2 + m03 t^3
        // y(t) = m10 + m11 t + m12 t^2 + m13 t^3
        //
        // First derivative (velocity vector):
        // x'(t) = m01 + 2 m02 t + 3 m03 t^2
        // y'(t) = m11 + 2 m12 t + 3 m13 t^2
        //
        // Second derivative (acceleration vector):
        // x''(t) = 2 m02 + 6 m03 t
        // y''(t) = 2 m12 + 6 m13 t
        //
        // Curvature.
        K0 = bc.curvature(0.0);
        plot_curvature = plot.create_connector(second_layer, line_style, plot_P0, plot_P0 + K0);

        double se = bc.stretching_energy(1e-9);
        double be = bc.bending_energy(1e-9);
        Dout(dc::notice, "stretching_energy = " << se);
        Dout(dc::notice, "bending_energy = " << be);
        Dout(dc::debug, "Total energy = " << (se + be));
      }
      auto curve = plot.create_curve(second_layer, curve_line_style, std::move(curve_points));

#if 0
      double radius = 1.0 / K0.norm();
      plot::Circle plot_curvature_circle;
      if (K0.norm() > 1e-9)
        plot_curvature_circle = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}),
            plot_P0 + radius * K0.direction(), radius);
#endif

      auto plot_extents = plot.create_rectangle(second_layer, rectangle_style, rect);
#endif

#if 1
      // Determine a quadratic Bezier going through P₀, P₁ and Pᵧ, tangent to D₀.
      //plot::Curve curve2;
      plot::BezierCurve plot_bezier_curve;
      plot::Connector plot_V_gamma;
      plot::Connector plot_S;
      plot::Connector plot_Q_gamma;
      plot::Point plot_P_gamma2;

      // BezierCurve works with math:: types.
      auto const& math_P0{plot_P0.raw()};
      auto const& math_P1{plot_P1.raw()};
      auto const& math_P_gamma{plot_P_gamma.raw()};
      auto const& math_D_gamma{D_gamma.raw()};

      do
      {
        BezierCurve qbc(math_P0, math_P1);
        double gamma;
        if (qbc.quadratic_from(math_P_gamma, math_D_gamma, point_gamma, gamma) && 0.0 < gamma && gamma < 1.0)
        {
          // Lets plot P_gamma using the gamma that was returned.
          Point const P_gamma{qbc.P(gamma)};
          plot_P_gamma2 = plot.create_point(second_layer, point_style({.color_index = 13, .filled_shape = 3}), P_gamma);

          Vector S = plot_P1 - plot_P_gamma;
          Vector Q_gamma = plot_P_gamma - plot_P0;
          Vector V_gamma = (1.0 - gamma) / gamma * Q_gamma + gamma / (1.0 - gamma) * S;
          plot_V_gamma = plot.create_connector(second_layer, solid_line_style, plot_P_gamma, plot_P_gamma + V_gamma);
          plot_S = plot.create_connector(second_layer, solid_line_style, plot_P_gamma, plot_P_gamma + S);
          plot_Q_gamma = plot.create_connector(second_layer, solid_line_style, plot_P_gamma, plot_P_gamma + Q_gamma);

#if 0
          std::vector<Point> curve_points2;
          for (int i = -400; i <= 800; ++i)
          {
            double t = i * 0.01;
            curve_points2.push_back(qbc.P(t));
          }
          curve2 = plot.create_curve(second_layer, curve_line_style, std::move(curve_points2));
#else
          plot_bezier_curve = plot.create_bezier_curve(second_layer, curve_line_style, qbc);
          Dout(dc::notice, "quadratic length = " << plot_bezier_curve.quadratic_arc_length());

#endif
        }
      }
      while (false);
#endif

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
