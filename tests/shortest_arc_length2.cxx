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

double angle_to_coordinate(double angle)
{
  return angle;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Shortest Arc Length Two Quadratic Beziers", 1200, 900);

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));

    auto reject_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("reject0_layer"));

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
        "Shortest arc length of two quadratic Beziers", {},
        "x", {},
        "y", {});
    plot.set_xrange({-M_PI, M_PI});
    plot.set_yrange({-M_PI, M_PI});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    int color_index2 = 3; //color_pool.get_and_use_color();
    Dout(dc::notice, "color_index2 = " << color_index2);
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    draw::LineStyle curve_line_style({.line_color = color::black, .line_width = 2.0});
    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});
    draw::LineStyle line_style({.line_color = color::gray, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::ArcStyle arc_style({.line_color = color::green, .line_width = 1.0});
    draw::ConnectorStyle connector_style(line_style({.line_color = color::coral, .dashes = {3.0, 3.0}}));

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, point_style, {-1.8849, -0.31415});
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, {0.6283, -0.31415});
    // Create a point P₂₂.
    auto plot_P2 = plot.create_point(second_layer, point_style, {1.811, 0.007943});

    // Create a point Q on the circle.
    double const circle_radius = 1.0;
    auto plot_Q0 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P0 + circle_radius * Direction{-0.628316});
    auto plot_Q1 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P1 + circle_radius * Direction{-2.08921});
    auto plot_Q2 = plot.create_point(second_layer, point_style({.color_index = 2}), plot_P2 + circle_radius * Direction{-0.933245});

    // Make all points draggable.
    window.register_draggable(plot, &plot_P0);
    window.register_draggable(plot, &plot_P1);
    window.register_draggable(plot, &plot_P2);

    window.register_draggable(plot, &plot_Q0, [&plot_P0, circle_radius](Point const& new_position)
        {
          return plot_P0 + circle_radius * Direction{plot_P0, new_position};
        }
    );
    window.register_draggable(plot, &plot_Q1, [&plot_P1, circle_radius](Point const& new_position)
        {
          return plot_P1 + circle_radius * Direction{plot_P1, new_position};
        }
    );
    window.register_draggable(plot, &plot_Q2, [&plot_P2, circle_radius](Point const& new_position)
        {
          return plot_P2 + circle_radius * Direction{plot_P2, new_position};
        }
    );

    std::vector<plot::Point> rejections;
    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      Dout(dc::notice, "P₀ = " << plot_P0);
      Dout(dc::notice, "P₁ = " << plot_P1);
      Dout(dc::notice, "P₂ = " << plot_P2);

      // Draw a circle around P₀.
      auto plot_circle0 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P0, circle_radius);
      // Draw a circle around P₁.
      auto plot_circle1 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P1, circle_radius);
      auto plot_circle2 = plot.create_circle(second_layer, solid_line_style({.line_color = color::gray}), plot_P2, circle_radius);

      // Draw a label for P₀, P₁ and P₂.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");
      auto P2_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P2, "P₂");

      // Draw line through P₀ and P₁.
      auto plot_line_P0P1 = plot.create_line(second_layer, line_style, plot::LineExtend::to, plot_P0, plot_P1);

      // Draw line through P₁ and P₂.
      auto plot_line_P1P2 = plot.create_line(second_layer, line_style, plot::LineExtend::to, plot_P1, plot_P2);

      Direction D0 = (plot_Q0 - plot_P0).direction();
      Direction D1 = (plot_Q1 - plot_P1).direction();
      Direction D2 = (plot_Q2 - plot_P2).direction();

      Dout(dc::notice, "D₀ = " << D0.as_angle());
      Dout(dc::notice, "D₁ = " << D1.as_angle());
      Dout(dc::notice, "D₂ = " << D2.as_angle());

      // Draw the tangent at P₀.
      auto plot_tangent_P0 = plot.create_line(second_layer, line_style, plot_P0, D0);

      // Draw the tangent in P₁.
      auto plot_tangent_P1 = plot.create_line(second_layer, line_style, plot_P1, D1);

      // Draw the tangent at P₂.
      auto plot_tangent_P2 = plot.create_line(second_layer, line_style, plot_P2, D2);

      // Calculate the angle offset at P₁.
      double angle_offset_1 = plot_line_P0P1.direction().as_angle() - plot_line_P1P2.direction().as_angle();
      Dout(dc::notice, "angle_offset_1 = " << angle_offset_1);

      // Draw the arc between P0P1 and the tangent in P₀.
      auto plot_arc01_0 = plot.create_arc(second_layer, arc_style, plot_P0, plot_line_P0P1.direction(), D0, 1.0);
      double arc01_0 = plot_arc01_0.end_angle() - plot_arc01_0.start_angle();

      // Draw the arc between P0P1 and the tangent in P₁.
      auto plot_arc01_1 = plot.create_arc(second_layer, arc_style, plot_P1, plot_line_P0P1.direction(), D1, 1.0);
      double arc01_1 = plot_arc01_1.end_angle() - plot_arc01_1.start_angle();

      // Draw the arc between P1P2 and the tangent in P₁.
      auto plot_arc12_1 = plot.create_arc(second_layer, arc_style({.line_color = color::red}), plot_P1, plot_line_P1P2.direction(), D1, 1.5);
      double arc12_1 = plot_arc12_1.end_angle() - plot_arc12_1.start_angle();

      // Draw the arc between P1P2 and the tangent in P₂.
      auto plot_arc12_2 = plot.create_arc(second_layer, arc_style({.line_color = color::red}), plot_P2, plot_line_P1P2.direction(), D2, 1.5);
      double arc12_2 = plot_arc12_2.end_angle() - plot_arc12_2.start_angle();

      Dout(dc::notice, "arc01_0 = " << arc01_0 << "; arc12_2 = " << arc12_2);
      Dout(dc::notice, "arc01_1 = " << arc01_1 << "; arc12_1 = " << arc12_1 << "; arc12_1 - arc01_1 = " << (arc12_1 - arc01_1));

      Vector Q1 = plot_P1 - plot_P0;
      Vector N1 = Q1.rotate_90_degrees();

      Dout(dc::notice, "Q1 = " << Q1);
      Dout(dc::notice, "N1 = " << N1);
      auto plot_arrow_Q1 = plot.create_connector(second_layer, connector_style({.line_color = color::blue}), plot_P0, plot_P0 + Q1);
      auto arrow_Q1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0 + 0.5 * Q1, "Q₁");
      auto plot_arrow_N1 = plot.create_connector(second_layer, connector_style({.line_color = color::blue}), plot_P0, plot_P0 + N1);
      auto arrow_N1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0 + 0.5 * N1, "N₁");

      bool generate_rejections = false; //rejections.empty();
      plot::BezierCurve plot_bezier_curve01;
      plot::Point plot_angles_point;

      if (generate_rejections)
        Debug(dc::notice.off());
      int const steps_per_pi = 250;
      int const steps = steps_per_pi;
      for (int a0 = -(steps - 1); a0 <= steps; ++a0)
      {
        double alpha0 = generate_rejections ? a0 * M_PI / steps_per_pi - plot_line_P0P1.direction().as_angle() : arc01_0;
        for (int a1 = -(steps - 1); a1 < steps; ++a1)
        {
          double alpha1 = generate_rejections ? a1 * M_PI / steps_per_pi - plot_line_P0P1.direction().as_angle() : arc01_1;

          double v0_div_q1 = 2.0 * std::sin(alpha1) / std::sin(M_PI + alpha0 - alpha1);
          Vector V0 = v0_div_q1 * std::cos(alpha0) * Q1 + v0_div_q1 * std::sin(alpha0) * N1;

          BezierCurve qbc01(plot_P0, plot_P1);
          qbc01.quadratic_from(V0);
#if 0
          double V0_dot_D0;
          double V1_dot_D1;
          if (generate_rejections)
          {
            Direction D0_(a0 * M_PI / steps_per_pi);
            Direction D1_(a1 * M_PI / steps_per_pi);
            V0_dot_D0 = V0.dot(D0_);
            V1_dot_D1 = qbc01.velocity(1.0).dot(D1_);
          }
          else
          {
            V0_dot_D0 = V0.dot(D0);
            V1_dot_D1 = qbc01.velocity(1.0).dot(D1);
          }
          bool reject = V0_dot_D0 <= 0.0 || V1_dot_D1 <= 0.0;
#else
          constexpr double to_loops = 1.0 / (2.0 * M_PI);
          double alpha0_in_loops = to_loops * alpha0;
          double alpha1_in_loops = to_loops * alpha1;
          double foo0 = alpha0_in_loops - std::floor(alpha0_in_loops);
          double foo1 = alpha1_in_loops - std::floor(alpha1_in_loops);
          bool reject = std::abs(foo0 - foo1) < 0.5;
#endif

          if (!generate_rejections)
          {
            plot_bezier_curve01 =
              plot.create_bezier_curve(second_layer, curve_line_style({.line_color = reject ? color::red : color::black}), qbc01);
            Dout(dc::notice, "Elastic energy 01 = " << utils::square(qbc01.quadratic_arc_length()));
            Dout(dc::notice, "Bending energy 01 = " << qbc01.bending_energy(1e-6));
            Dout(dc::notice, "Quadratic E    01 = " << qbc01.quadratic_bending_energy());
          }

          // Plot rejection.
          double x = angle_to_coordinate(alpha0);
          double y = angle_to_coordinate(alpha1);

          if (generate_rejections)
          {
            Dout(dc::notice, "Adding point with reject=" << reject << " at " << x << ", " << y);
            rejections.push_back(plot.create_point(reject_layer, point_style({.color_index = reject ? 13 : 10, .filled_shape = 3}), {x, y}));
          }
          else
          {
            plot_angles_point = plot.create_point(reject_layer, point_style({.color_index = 5, .filled_shape = 3}), {x, y});
            break;
          }
        }
        if (!generate_rejections)
          break;
      }
      if (generate_rejections)
        Debug(dc::notice.on());

      plot::BezierCurve plot_bezier_curve12;
      BezierCurve qbc12(plot_P1, plot_P2);
      if (qbc12.quadratic_from(D1, D2))
      {
        plot_bezier_curve12 = plot.create_bezier_curve(second_layer, curve_line_style, qbc12);
        Dout(dc::notice, "Elastic energy 12 = " << utils::square(qbc12.quadratic_arc_length()));
        Dout(dc::notice, "Bending energy 12 = " << qbc12.bending_energy(1e-6));
        Dout(dc::notice, "Quadratic E    12 = " << qbc12.quadratic_bending_energy());
      }

      // Draw V0.
//      auto plot_V0 = plot.create_connector(second_layer, connector_style, plot_P0, plot_P0 + 0.5 * qbc01.V0());

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // Block until the user moved a draggable object, then go to the top of loop for a redraw.
      window.handle_dragging();
    }

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
