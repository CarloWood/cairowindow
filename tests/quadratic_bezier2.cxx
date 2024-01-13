#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Vector.h"
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
        "Relationship between P₀, P₁, θ, w, v and s", {},
        "x", {},
        "y", {});
    plot.set_xrange({-10, 11});
    plot.set_yrange({-10.5, 10.5});
    plot.add_to(background_layer, true);

    utils::ColorPool<32> color_pool;
    int color_index = color_pool.get_and_use_color();
    int color_index2 = 3; //color_pool.get_and_use_color();
    Dout(dc::notice, "color_index2 = " << color_index2);
    int filled_shape = 1;
    draw::PointStyle point_style(color_index, filled_shape);
    draw::TextStyle<> label_style{.position = draw::centered_left_of, .font_size = 18.0, .offset = 10};
    draw::LineStyle curve_line_style{.line_width = 1.0};
    draw::LineStyle solid_line_style{.line_color = color::black, .line_width = 1.0};
    draw::LineStyle line_style{.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}};
    draw::ArcStyle arc_style{.line_color = color::blue, .line_width = 1.0};

    double a = -1.0;
    for (int j = 0;; ++j)
    {
      double theta = j * M_PI / 32;

      // P₀, the point at t=0.
      Point P0(0.0, 0.0);
      auto plot_P0 = plot.create_point(second_layer, P0, point_style);
      auto P0_label = plot.create_text(second_layer, P0, "P₀", label_style);

      // P₁, the point at t=1.
      Point P1(1.0, 0.0);
      auto plot_P1 = plot.create_point(second_layer, P1, point_style);
      auto P1_label = plot.create_text(second_layer, P1, "P₁", label_style({.position = draw::centered_right_of}));

      // Draw a cirle around the midpoint of P₀P₁ with radius |P₀P₁|/2.
      Vector P0P1(P0, P1);
      Point P0P1_circle_center = P0 + 0.5 * P0P1;
      auto plot_P0P1_circle_center = plot.create_point(second_layer, P0P1_circle_center, point_style);
      double P0P1_circle_radius = 0.5 * P0P1.length();
      auto plot_P0P1_circle = plot.create_circle(second_layer, P0P1_circle_center, P0P1_circle_radius,
          line_style({.line_color = color::gray}));

      // Define the matrix M.
#if 0
      double m00 = std::cos(theta);
      double m01 = 0.0;
      double m10 = -std::sin(theta) - a;
      double m11 = a;
#endif
      double m00 = 1.0 + a * std::sin(theta);
      double m01 = -a * std::sin(theta);
      double m10 = -a * std::cos(theta);
      double m11 = a * std::cos(theta);

      auto xt = [=](double t){ return P0.x() + t * (m00 + m01 * t); };
      auto yt = [=](double t){ return P0.y() + t * (m10 + m11 * t); };

      // Let t run from v-4 to v+4, and then plot P₀ + t M [1 t].
      double v = 0; //0.5 * (1.0 + std::sin(theta) / a);
      std::vector<Point> curve_points;
      for (double t = v - 40.0; t <= v + 40.0; t += 0.1)
      {
        double x = xt(t);
        double y = yt(t);
        curve_points.emplace_back(x, y);
      }
      auto curve = plot.create_curve(second_layer, std::move(curve_points), curve_line_style);

      // Draw a line through P₀ and P₁.
      auto line_through_P0_and_P1 = plot.create_line(second_layer, P0, P1, solid_line_style);

      window.set_send_expose_events(true);

      std::this_thread::sleep_for(std::chrono::milliseconds(20));
      std::cin.get();

      window.set_send_expose_events(false);
    }

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
