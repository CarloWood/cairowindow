#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Plot.h"
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
    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "Lachance-Schwartz", {},
        "x", {},
        "y", {});
    plot.set_xrange({-10, 10});
    plot.set_yrange({-10, 10});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style(color_pool.get_and_use_color(), 1);
    draw::TextStyle<> label_style{.position = draw::centered_left_of, .font_size = 18.0, .offset = 10};
    draw::LineStyle curve_line_style{.line_width = 1.0};

    // Create a point Q₁.
    auto Q1 = plot.create_point(second_layer, {-2.0, -1.0}, point_style);
    // Create a point Q₂.
    auto Q2 = plot.create_point(second_layer, {2.0, -1.0}, point_style);
    // Create a point Q₃.
    auto Q3 = plot.create_point(second_layer, {0.0, 1.0}, point_style);
    // Create a point Q₄.
    auto Q4 = plot.create_point(second_layer, {0.0, -2.0}, point_style);

    // Make all points draggable.
    window.register_draggable(plot, &Q1);
    window.register_draggable(plot, &Q2);
    window.register_draggable(plot, &Q3);
    window.register_draggable(plot, &Q4);

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a label for Q₁, Q₂, Q₃ and Q₄.
      auto Q1_label = plot.create_text(second_layer, Q1, "Q₁", label_style);
      auto Q2_label = plot.create_text(second_layer, Q2, "Q₂", label_style);
      auto Q3_label = plot.create_text(second_layer, Q3, "Q₃", label_style);
      auto Q4_label = plot.create_text(second_layer, Q4, "Q₄", label_style);

      Eigen::Matrix3d R;
      R << Q1.x(), Q1.y(), 1.0,
           Q2.x(), Q2.y(), 1.0,
           Q3.x(), Q3.y(), 1.0;

      Eigen::RowVector3d Q4_1;
      Q4_1 << Q4.x(), Q4.y(), 1.0;

      auto q = Q4_1 * R.inverse();

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
        window.handle_dragging();
        continue;
      }

      std::array<plot::Curve, 2> curve;
      for (int solution = 0; solution < 2; ++solution)
      {
        double discriminant = utils::square(-2.0 * q(1) * q(2)) - 4.0 * q(1) * (1.0 - q(1)) * q(2) * (1.0 - q(2));
        double a = (2.0 * q(1) * q(2) + ((solution == 0) ? -1 : 1) * std::sqrt(discriminant)) / (2.0 * q(1) * (1.0 - q(1)));

        // Only draw parabola's that have Q2 before Q0 and Q1.
        if (a > 0.0)
          continue;

        Eigen::Matrix3d V;
        V << 0.0, 0.0, 1.0,
             a*a,   a, 1.0,
             1.0, 1.0, 1.0;

        Eigen::Matrix3d VinvR = V.inverse() * R;

        std::vector<Point> curve_points;
        {
          for (int i = -200; i <= 200; ++i)
          {
            double t = i * 0.1;
            Eigen::RowVector3d T;
            T << (t * t), t, 1.0;
            auto Phi = T * VinvR;
            curve_points.emplace_back(Phi(0), Phi(1));
          }
        }
        curve[solution] = plot.create_curve(second_layer, std::move(curve_points),
            curve_line_style({.line_color = (solution == 0) ? color::green : color::blue}));
      }

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
