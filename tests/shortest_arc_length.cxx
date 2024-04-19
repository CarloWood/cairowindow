#include "sys.h"
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
    Window window("Shortest Arc Length Quadratic Bezier", 1200, 900);

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
        "Shortest arc length of quadratic Bezier", {},
        "x", {},
        "y", {});
    plot.set_xrange({-10, 10});
    plot.set_yrange({-10, 10});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    int color_index2 = 3; //color_pool.get_and_use_color();
    Dout(dc::notice, "color_index2 = " << color_index2);
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});
    draw::LineStyle line_style({.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::ArcStyle arc_style({.line_color = color::blue, .line_width = 1.0});
    draw::ConnectorStyle connector_style(line_style({.line_color = color::coral, .dashes = {3.0, 3.0}}));

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, point_style, {-6.0, -1.0});
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, {2.0, 0.0});
    // Create a point Pᵧ.
    auto plot_P_gamma = plot.create_point(second_layer, point_style({.filled_shape = 10}), {-4.0, 1.0});

    // Make all points draggable.
    window.register_draggable(plot, &plot_P0);
    window.register_draggable(plot, &plot_P1);
    window.register_draggable(plot, &plot_P_gamma);

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a label for P₀, P₁ and Pᵧ.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");
      auto P_gamma_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P_gamma, "Pᵧ");

      // Draw line through P₀ and Pᵧ.
      auto plot_line_P0Pgamma = plot.create_line(second_layer, line_style, plot::LineExtend::both, plot_P0, plot_P_gamma);

      // Draw line through P₁ and Pᵧ.
      auto plot_line_P1Pgamma = plot.create_line(second_layer, line_style, plot::LineExtend::both, plot_P_gamma, plot_P1);

      Direction const D_gamma((Vector{plot_line_P0Pgamma.direction()} + Vector{plot_line_P1Pgamma.direction()}).direction());

      plot::BezierCurve plot_bezier_curve;
      plot::Connector plot_velocity_gamma;
      plot::Line plot_tangent_P0;
      plot::Line plot_tangent_P1;
      plot::Line plot_normal_P_gamma;
      plot::Point plot_C;
      do
      {
        BezierCurve qbc(plot_P0, plot_P1);
        double gamma;
        if (qbc.quadratic_from(plot_P_gamma, D_gamma, point_gamma, gamma) && 0.0 < gamma && gamma < 1.0)
        {
          plot_bezier_curve = plot.create_bezier_curve(second_layer, curve_line_style, qbc);
          Dout(dc::notice, "quadratic length = " << qbc.quadratic_arc_length());

          // Draw the velocity vector at Pᵧ.
          plot_velocity_gamma = plot.create_connector(second_layer,
              solid_line_style({.line_color = color::green}),
              plot_P_gamma, plot_P_gamma + plot_bezier_curve.velocity(gamma));

          // Draw a line through P₀ tangent to the curve.
          plot_tangent_P0 = plot.create_line(second_layer, line_style({.line_color = color::cyan}), plot_P0, qbc.V0().direction());
          // Draw a line through P₁ tangent to the curve.
          plot_tangent_P1 = plot.create_line(second_layer, line_style({.line_color = color::cyan}), plot_P1, qbc.velocity(1.0).direction());
          // Draw a point at their intersection.
          Point C = plot_tangent_P0.intersection_with(plot_tangent_P1);
          plot_C = plot.create_point(second_layer, point_style({.color_index = 4}), C);
          // Draw a line through Pᵧ, perpendicular to the curve.
          plot_normal_P_gamma = plot.create_line(second_layer, solid_line_style, plot_P_gamma, plot_velocity_gamma.direction().normal());
        }
      }
      while (false);

      Dout(dc::notice, "cubic length = " << plot_bezier_curve.arc_length(0.001));
      Dout(dc::notice, "quadratic length = " << plot_bezier_curve.quadratic_arc_length());

      // Draw V, the parabola vertex point, which resides at t=v.
      BezierCurveMatrix const& m = plot_bezier_curve.M();
      double v = -0.5 * m.coefficient[1].dot(m.coefficient[2]) / m.coefficient[2].length_squared();
      plot::Point V;
      plot::Text V_label({0, 0}, "");
      if (!std::isnan(v))
      {
        V = plot.create_point(second_layer, point_style, plot_bezier_curve.P(v));
        V_label = plot.create_text(second_layer, label_style({.position = draw::centered_below}), V, "V");
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
