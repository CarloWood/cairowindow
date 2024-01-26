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
    Window window("Quadratic Bezier", 1200, 900);

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

    utils::ColorPool<32> color_pool;
    int color_index = color_pool.get_and_use_color();
    int color_index2 = 3; //color_pool.get_and_use_color();
    Dout(dc::notice, "color_index2 = " << color_index2);
    int filled_shape = 1;
    draw::PointStyle point_style(color_index, filled_shape);
    draw::TextStyle<> label_style{.position = draw::centered_left_of, .font_size = 18.0, .offset = 10};
    draw::TextStyle<> slider_style{.position = draw::centered_below, .font_size = 18.0, .offset = 10};
    draw::LineStyle curve_line_style{.line_width = 1.0};
    draw::LineStyle solid_line_style{.line_color = color::black, .line_width = 1.0};
    draw::LineStyle line_style{.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}};
    draw::ArcStyle arc_style{.line_color = color::blue, .line_width = 1.0};

    // Draw a slider for alpha.
    auto slider_alpha = plot.create_slider(second_layer, {928, 83, 7, 400}, 1.0, 0.1, 10.0);
    auto slider_alpha_label = plot.create_text(second_layer, Pixel{928, 483}, "α", slider_style);

    // Draw a slider for m03.
    auto slider_m03 = plot.create_slider(second_layer, {978, 83, 7, 400}, 0.0, -1.0, 1.0);
    auto slider_m03_label = plot.create_text(second_layer, Pixel{978, 483}, "m03", slider_style);

    // Draw a slider for m13.
    auto slider_m13 = plot.create_slider(second_layer, {1028, 83, 7, 400}, 0.0, -1.0, 3.0);
    auto slider_m13_label = plot.create_text(second_layer, Pixel{1028, 483}, "m13", slider_style);

    // Create a point P₀.
    auto plot_P0 = plot.create_point(second_layer, {-5.0, -3.0}, point_style);
    // Create a point P₁.
    auto plot_P1 = plot.create_point(second_layer, {3.0, -1.0}, point_style);

    // Create a point Q on a circle around P₀.
    double P0_circle_radius = 2.0;
    double phi = M_PI / 10.0;
    Direction V0_dir{phi};
    auto plot_Q = plot.create_point(second_layer, plot_P0 + P0_circle_radius * V0_dir, point_style);

    // Make P₀ and P₁ draggable.
    window.register_draggable(plot, &plot_P0);
    window.register_draggable(plot, &plot_P1);
    // Make Q draggable along the circle.
    window.register_draggable(plot, &plot_Q, [&plot_P0, P0_circle_radius](Point const& new_position)
        {
          return plot_P0 + P0_circle_radius * Direction{plot_P0, new_position};
        }
    );

    Point prev_Q = plot_Q;

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Draw a label for P₀ and P₁.
      auto P0_label = plot.create_text(second_layer, plot_P0, "P₀", label_style({.position = draw::centered_right_of}));
      auto P1_label = plot.create_text(second_layer, plot_P1, "P₁", label_style({.position = draw::centered_right_of}));

      // Draw line through P₀ and P₁.
      auto plot_line_P0P1 = plot.create_line(second_layer, plot_P0, plot_P1, line_style, plot::LineExtend::both);

      // Draw a draggable point (Q) on a circle around P₀.
      auto plot_P0_circle = plot.create_circle(second_layer, plot_P0, P0_circle_radius,
          line_style({.line_color = color::gray}));

      Dout(dc::notice, "V0_dir = " << V0_dir);
      if (prev_Q != plot_Q)
      {
        Dout(dc::notice, "prev_Q != plot_Q");
        prev_Q = plot_Q;
        V0_dir = (plot_Q - plot_P0).direction();
      }
      else
      {
        Dout(dc::notice, "prev_Q == plot_Q");
        // Update Q in case P₀ was moved.
        plot_Q = plot.create_point(second_layer, plot_P0 + P0_circle_radius * V0_dir, point_style);
        // Make Q draggable along the circle.
        window.register_draggable(plot, &plot_Q, [&plot_P0, P0_circle_radius](Point const& new_position)
            {
              return plot_P0 + P0_circle_radius * Direction{plot_P0, new_position};
            }
        );
        prev_Q = plot_Q;
      }
      Dout(dc::notice, "V0_dir = " << V0_dir);

      // Draw an arrow from P₀ to Q.
      auto plot_Q_arrow = plot.create_connector(second_layer, plot_P0, plot_Q, line_style({.line_color = color::orange}));

      // Define the matrix.
      //      ⎡Px₀  α·Vx₀  m₀₂  Px₁-(Px₀+α·Vx₀+m₀₂)⎤
      //  M = ⎣Py₀  α·Vy₀  m₁₂  Py₁-(Py₀+α·Vy₀+m₁₂)⎦
      double alpha = slider_alpha.value();
      Vector V0 = plot_Q_arrow;
      double m00 = plot_P0.x();
      double m10 = plot_P0.y();
      double m01 = alpha * V0.x();
      double m11 = alpha * V0.y();
      double m03 = slider_m03.value();
      double m13 = slider_m13.value();
      double m02 = plot_P1.x() - (m00 + m01 + m03);
      double m12 = plot_P1.y() - (m10 + m11 + m13);

      auto xt = [=](double t){ return m00 + t * (m01 + t * (m02 + t * m03)); };
      auto yt = [=](double t){ return m10 + t * (m11 + t * (m12 + t * m13)); };

      std::vector<Point> curve_points;
      for (int i = -100; i <= 100; ++i)
      {
        double t = i * 0.1;
        curve_points.emplace_back(xt(t), yt(t));
      }
      auto curve = plot.create_curve(second_layer, std::move(curve_points), curve_line_style);

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
