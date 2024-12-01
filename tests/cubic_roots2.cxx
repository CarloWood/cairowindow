#include "sys.h"
#include "math/CubicPolynomial.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include "utils/square.h"
#include <thread>
#include <cmath>
#include "debug.h"

bool almost_equal(double x1, double x2)
{
  return std::abs(x1 - x2) < 1e-6;
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
    Window window("Cubic roots", 1200, 1000);

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
        "Cubic roots", {},
        "x", {},
        "y", {});
    plot.set_xrange({-10.0, 10.0});
    plot.set_yrange({-100.0, 100.0});
    plot.add_to(background_layer, false);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});
    draw::LineStyle dashed_line_style({.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    draw::PointStyle root0_style({.color_index = 12, .filled_shape = 3});
    draw::PointStyle root1_style({.color_index = 2, .filled_shape = 3});
    draw::PointStyle root2_style({.color_index = 10, .filled_shape = 3});
    draw::PointStyle inflection_style({.color_index = 13, .filled_shape = 5});
    draw::CircleStyle circle_style({.position = draw::at_center});

    // Draw a slider for r₀.
    auto slider_r0 = plot.create_slider(second_layer, {1028, 83, 7, 400}, 0.0001, -10.0, 10.0);
    auto slider_r0_label = plot.create_text(second_layer, slider_style, Pixel{1028, 483}, "r₀");

    // Draw a slider for r₁.
    auto slider_r1 = plot.create_slider(second_layer, {1078, 83, 7, 400}, 1.5, -10.0, 10.0);
    auto slider_r1_label = plot.create_text(second_layer, slider_style, Pixel{1078, 483}, "r₁");

    // Draw a slider for r₂.
    auto slider_r2 = plot.create_slider(second_layer, {1128, 83, 7, 400}, 2.5, -10.0, 10.0);
    auto slider_r2_label = plot.create_text(second_layer, slider_style, Pixel{1128, 483}, "r₂");

    // Draw x and y axis.
    auto x_axis = plot.create_line(background_layer, solid_line_style({.line_color = color::slategrey}), Point{0.0, 0.0}, Direction::left);
    auto y_axis = plot.create_line(background_layer, solid_line_style({.line_color = color::slategrey}), Point{0.0, 0.0}, Direction::up);

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Extract slider values.
      double r0 = slider_r0.value();
      double r1 = slider_r1.value();
      double r2 = slider_r2.value();

      // Correct slider values.
      if (std::abs(r0) > std::abs(r1))
      {
        r1 = std::copysign(std::abs(r0), r1);
        slider_r1.set_value(r1);
      }
      if (std::abs(r1) > std::abs(r2))
      {
        r2 = std::copysign(std::abs(r1), r2);
        slider_r2.set_value(r2);
      }

      // Calculate the cubics coefficients.
      double c0 = -r0 * r1 * r2;
      double c1 = r0 * r1 + r0 * r2 + r1 * r2;
      double c2 = -(r0 + r1 + r2);
      math::CubicPolynomial cubic(c0, c1, c2, 1.0);
      Dout(dc::notice, cubic);

      double const d = std::max(0.0, utils::square(c2) - 3.0 * c1);     // Prevent floating-point round off errors to result in a negative d.
      double x_max = (-c2 - std::sqrt(d)) / 3.0;
      double x_min = (-c2 + std::sqrt(d)) / 3.0;
      ASSERT(x_max <= x_min);
      double P0 = cubic(0.0);
      double Pxmax = cubic(x_max);
#if 0
      if (x_max > 0)
      {
        if (2 * Pxmax > (3 * std::abs(P0) + P0))
          Dout(dc::notice, "Class A");
        else
          Dout(dc::notice, "Class B");
      }
      else
      {
        if (2 * Pxmax > (3 * std::abs(P0) + P0))
          Dout(dc::notice, "Class D");
        else
          Dout(dc::notice, "Class C");
      }
#endif

      // Draw the cubic.
      plot::BezierFitter plot_cubic;
      plot_cubic.solve(
          [&](double x) -> Point { return {x, cubic(x)}; },
          plot.viewport());
      plot.add_bezier_fitter(second_layer, solid_line_style, plot_cubic);

      // Plot the roots, using Red for r₀, Green for r₁ and Blue for r₂.
      auto plot_r0 = plot.create_point(second_layer, root0_style, {r0, 0.0});
      auto plot_r1 = plot.create_point(second_layer, root1_style, {r1, 0.0});
      auto plot_r2 = plot.create_point(second_layer, root2_style, {r2, 0.0});

      // Draw a vertical line at x_max and x_min.
      auto plot_x_max = plot.create_line(second_layer, dashed_line_style({.line_color = x_max < 0.0 ? color::brown : color::lime}),
          Point{x_max, 0.0}, Direction::up);
      auto plot_x_min = plot.create_line(second_layer, dashed_line_style({.line_color = x_min < 0.0 ? color::brown : color::lime}),
          Point{x_min, 0.0}, Direction::up);
      // Draw a vertical line at the inflection point.
      double x_inflection_point = -c2 / 3.0;
      auto plot_inflection = plot.create_line(second_layer, dashed_line_style, Point{x_inflection_point, 0.0}, Direction::up);
      // Evaluate the cubic at the inflection point.
      double y_inflection_point = cubic(x_inflection_point);
      auto plot_inflection_point = plot.create_point(second_layer, inflection_style, {x_inflection_point, y_inflection_point});

      Dout(dc::notice, "x_max = " << x_max << ", x_min = " << x_min << ", r0 = " << r0 << ", r1 = " << r1 << ", r2 = " << r2);

      // Sort the known roots to what we can choose an algorithm for:
      // * left: left than x_max.
      // * middle: in between x_max and x_min.
      // * right: right of x_min.
      std::array<double, 3> roots{{r0, r1, r2}};
      std::sort(roots.begin(), roots.end());
      constexpr int left = 0;
      constexpr int middle = 1;
      constexpr int right = 2;

      // Which root would we calculate?
      // If inflection point is in the upper-right or lower-left quadrant
      // we'll pick the same root as if 0 < x_max < x_min or x_max < x_min < 0 respectively.
      //  2 | 0
      //  --+--
      //  3 | 1
      int inflection_quadrant = (x_inflection_point < 0.0 ? 2 : 0) | (y_inflection_point < 0.0 ? 1 : 0);
      bool changed = false;
      double root;
      if (0.0 < x_max)
      {
        // Calculate the root that is less than x_max.
        root = roots[left];
      }
      else if (x_min < 0.0)
      {
        // Calculate the root that is larger than x_min.
        root = roots[right];
      }
      else if (inflection_quadrant == 0 || inflection_quadrant == 3)
      {
        changed = true;
        // Pick the root that is on the outer side of the local extreme that is closest to 0.
        if (std::abs(x_max) < std::abs(x_min) && std::abs(x_max) < std::abs(x_inflection_point))
          root = roots[left];
        else if (std::abs(x_min) < std::abs(x_inflection_point))
          root = roots[right];
        else
          root = roots[middle];
      }
      else
      {
        ASSERT(x_max <= 0.0 && 0.0 <= x_min);
        // Calculate the root that is larger than x_max but less than x_min.
        root = roots[middle];
      }
      auto plot_root = plot.create_connector(second_layer,
          solid_line_style({.line_color = changed ? color::red : color::black, .line_width = 2.0}),
          Point{root, cubic(root) + 10.0}, Point{root, cubic(root) + 2.0});

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
