#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Plot.h"
#include "cairowindow/BezierCurve.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include "utils/almost_equal.h"
#include <Eigen/Dense>
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
    Window window("Four point parabolic interpolation", 1200, 900);

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
        "Lachance-Schwartz", {},
        "x", {},
        "y", {});
    plot.set_xrange({-2, 6});
    plot.set_yrange({-4, 4});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    draw::PointStyle point_circle_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 10});
    draw::PointStyle point_square_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 0});
    draw::PointStyle point_triangle_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 7});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::LineStyle line_style({.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::LineStyle curve_line_style({.line_width = 1.0});

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, point_style, {-0.3, -0.6});
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, point_style, {1.0, 0.4});
    // Create a point Pᵦ.
    auto plot_P_beta = plot.create_point(second_layer, point_style, {-1.0, 2.0});
    // Create a point Pᵧ.
    auto plot_P_gamma = plot.create_point(second_layer, point_style, {2.2, 3.3});

    // BezierCurve uses math:: types.
    auto const& math_P0{plot_P0.raw()};
    auto const& math_P1{plot_P1.raw()};
    auto const& math_P_beta{plot_P_beta.raw()};
    auto const& math_P_gamma{plot_P_gamma.raw()};

    // Make all points draggable.
    window.register_draggable(plot, &plot_P0);
    window.register_draggable(plot, &plot_P1);
    window.register_draggable(plot, &plot_P_beta);
    window.register_draggable(plot, &plot_P_gamma);

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a label for P₀, P₁, Pᵦ and Pᵧ.
      auto P0_label = plot.create_text(second_layer, label_style, plot_P0, "P₀");
      auto P1_label = plot.create_text(second_layer, label_style, plot_P1, "P₁");
      auto P_beta_label = plot.create_text(second_layer, label_style, plot_P_beta, "Pᵦ");
      auto P_gamma_label = plot.create_text(second_layer, label_style, plot_P_gamma, "Pᵧ");

      Eigen::Matrix3d R;
      R << plot_P0.x(), plot_P0.y(), 1.0,
           plot_P_beta.x(), plot_P_beta.y(), 1.0,
           plot_P1.x(), plot_P1.y(), 1.0;

      Eigen::RowVector3d P_gamma_1;
      P_gamma_1 << plot_P_gamma.x(), plot_P_gamma.y(), 1.0;

      auto q = P_gamma_1 * R.inverse();

      // See https://deepblue.lib.umich.edu/bitstream/handle/2027.42/29347/0000415.pdf
      int count = 0;
      for (int i = 0; i < 3; ++i)
      {
        if (utils::almost_equal(q(i), 1.0, 10e-9))
          ++count;
      }
      if (count >= 2 || q(0) * q(1) * q(2) >= 0.0)
      {
        // There is no parabola going through all four points.
        window.set_send_expose_events(true);
        // Block till the next user interaction.
        if (window.handle_input_events())
          continue;
        // The program was terminated by the user.
        break;
      }

      Dout(dc::notice, "q = " << q);

      std::array<plot::Point, 2> marker;
      std::array<plot::Connector, 2> control;
      plot::Curve curve2;

      BezierCurve bc(math_P0, math_P1);
      if (bc.quadratic_from(math_P_beta, math_P_gamma))
      {
        std::vector<Point> curve_points;
        {
          for (int i = -200; i <= 300; ++i)
          {
            double t = i * 0.01;
            Point Pt{bc.P(t)};
            curve_points.push_back(Pt);
          }
        }
        curve2 = plot.create_curve(second_layer, curve_line_style({.line_color = color::red}), std::move(curve_points));

        Point const P0{bc.P0()};
        Point const P1{bc.P1()};
        marker[0] = plot.create_point(second_layer, point_circle_style, P0);
        marker[1] = plot.create_point(second_layer, point_square_style, P1);

        Point const C0{bc.C0()};
        Point const C1{bc.C1()};
        control[0] = plot.create_connector(second_layer, line_style, plot_P0, C0);
        control[1] = plot.create_connector(second_layer, line_style, plot_P1, C1);
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
