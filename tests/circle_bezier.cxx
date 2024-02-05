#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Plot.h"
#include "cairowindow/BezierFitter.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
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
    Window window("Circle Bezier", 1200, 900);

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
        "Circle Bezier test", {},
        "x", {},
        "y", {});
    plot.set_xrange({0, 200});
    plot.set_yrange({0, 200});
    plot.add_to(background_layer, true);

    draw::PointStyle point_style({.color_index = 0, .filled_shape = 1});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::LineStyle line_style({.line_color = color::red, .line_width = 1.0});
    draw::BezierCurveStyle bezier_curve_style({.line_color = Color::next_color(), .line_width = 1.0});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});

    auto slider_velocity = plot.create_slider(second_layer, {978, 83, 7, 400}, 105.0, 0.0, 150.0);
    auto slider_velocity_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "v");

    auto plot_circle = plot.create_circle(background_layer, line_style, Point{100.0, 100.0}, 80.0);

#if 1
    //BezierFitter fitter([](double t) -> Point { return {100.0 + 80.0 * std::cos(t), 100.0 + 80.0 * std::sin(t)}; }, {-M_PI, M_PI}, {0.0, 0.0, 200.0, 200.0}, 0.1);

    BezierFitter fitter([](double t) -> Point {
        return {100.0 * t / M_PI, 100.0 + 100.0 * std::sin(5.0 * t * t)};
    }, {0, 2.0 * M_PI}, {0.0, 0.0, 200.0, 200.0}, 0.001);

    std::vector<BezierCurve> result = fitter.solve();
    std::vector<plot::BezierCurve> curves(result.size());
    std::vector<plot::Point> points0(result.size());
    int i = 0;
    for (auto&& bezier : result)
    {
      points0[i] = plot.create_point(second_layer, point_style, result[i].P(0));
      curves[i] = plot.create_bezier_curve(second_layer, bezier_curve_style, bezier);
      ++i;
    }
#else
    Point P0{20.0, 100};
    Point P1{180.0, 100};

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      double v = slider_velocity.value();
      Point C0 = P0 + Vector{0.0, v};
      Point C1 = P1 + Vector{0.0, v};
      BezierCurve bezier(P0, C0, C1, P1);

#if 0
      std::vector<Point> curve_points;
      for (int i = 0; i <= 100; ++i)
      {
        double t = i * 0.01;
        curve_points.push_back(bezier.P(t));
      }
      auto plot_curve = plot.create_curve(second_layer, curve_line_style, std::move(curve_points));
#else
      auto plot_curve = plot.create_bezier_curve(second_layer, bezier_curve_style, bezier);
#endif

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // Block until the user moved a draggable object, then go to the top of loop for a redraw.
      window.handle_dragging();
    }
#endif

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
