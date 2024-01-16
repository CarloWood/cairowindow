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
        "Relationship between P₀, P₁, θ, w, v and s", {},
        "x", {},
        "y", {});
    plot.set_xrange({-1, 2});
    plot.set_yrange({-1.5, 1.5});
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

    double rsw = 0.5;           // This is probably the reciprocal of the square of the "width" (see w in README.bezier).
    for (int j = 0;; ++j)
    {
      double theta = j * M_PI / 32;

      // Create rotation matrix for rotating -theta.
      Matrix R{std::cos(theta), std::sin(theta), -std::sin(theta), std::cos(theta)};

      // P₀, the point at t=0, was translated to the origin.
      auto plot_P0 = plot.create_point(second_layer, {0, 0}, point_style);
      auto P0_label = plot.create_text(second_layer, plot_P0, "P₀", label_style);

      // P₁, the point at t=1.
      Point P1(1.0, 0.0);
      auto plot_P1 = plot.create_point(second_layer, P1, point_style);
      auto P1_label = plot.create_text(second_layer, P1, "P₁", label_style({.position = draw::centered_right_of}));

      // Draw a line through P₀ and P₁.
      auto line_through_P0_and_P1 = plot.create_line(second_layer, plot_P0, P1, solid_line_style);

      // Draw a cirle around the midpoint of P₀P₁ with radius |P₀P₁|/2.
      Vector P0P1(P1);
      Point P0P1_circle_center = (0.5 * P0P1).point();
      auto plot_P0P1_circle_center = plot.create_point(second_layer, P0P1_circle_center, point_style);
      double P0P1_circle_radius = 0.5 * P0P1.length();
      auto plot_P0P1_circle = plot.create_circle(second_layer, P0P1_circle_center, P0P1_circle_radius,
          line_style({.line_color = color::gray}));

      // Pᵦ, a point at t < 0.
      double x_beta = -0.5;
      double y_beta = 0.25;
      Point P_beta(x_beta, y_beta);
      auto plot_P_beta = plot.create_point(second_layer, P_beta, point_style);
      auto P_beta_label = plot.create_text(second_layer, P_beta, "Pᵦ", label_style({.position = draw::centered_right_of}));

      // Pᵧ, a point at t > 1.
      double x_gamma = 1.8;
      double y_gamma = 0.15;
      Point P_gamma(x_gamma, y_gamma);
      auto plot_P_gamma = plot.create_point(second_layer, P_gamma, point_style);
      auto P_gamma_label = plot.create_text(second_layer, P_gamma, "Pᵧ", label_style({.position = draw::centered_right_of}));

      double y_span = y_gamma - y_beta;
#if 0
      // Helper variables.
      double x_span = x_gamma - x_beta;
      double subexpr50 = 4.0 * y_gamma * y_span * (x_gamma * (1.0 - x_gamma) - y_gamma * (1.0 - x_beta) * x_beta / y_beta);
      double tan_theta = std::sqrt(utils::square(2.0 * x_span * y_gamma) + subexpr50) - x_span / y_span;
#endif

      double subexpr = 0.5 * y_beta * (-2.0 * (x_gamma - x_beta) * y_gamma +
          std::sqrt(utils::square(-2.0 * x_beta + 2.0 * x_gamma) * utils::square(y_gamma) -
            4.0 * y_gamma * y_span * (-x_gamma + utils::square(x_gamma) +
              (x_beta * y_gamma)/y_beta - (utils::square(x_beta) * y_gamma)/y_beta)));

      double z = x_beta + subexpr / (y_gamma * y_span);
      // The t^2 factor of x(t):
      double m01 = subexpr / (y_gamma * y_span * (1.0 - z) * z);

//      double z = x_beta + y_beta * tan_theta;
      // Define the matrix M.
      double m10 = y_beta / (z * (1.0 - z));
//      double m01 = m10 * tan_theta;
      double m11 = -m10;
      double m00 = 1.0 - m01;

#if 0   // The three point case - with free theta.
      // Helper variable.
      double z = x_beta + y_beta * std::tan(theta);
      // Define the matrix M.
      double m11 = y_beta / (z * (z - 1));
      double m10 = -m11;
      double m01 = m10 * std::tan(theta);
      double m00 = 1.0 - m01;
#endif

      auto xt = [=](double t){ return t * (m00 + m01 * t); };
      auto yt = [=](double t){ return t * (m10 + m11 * t); };

      std::vector<Point> curve_points;
      for (int i = -100; i <= 40; ++i)
      {
        double t = i * 0.1;
        Vector v{xt(t), yt(t)};
        curve_points.push_back(v.point());
      }
      auto curve = plot.create_curve(second_layer, std::move(curve_points), curve_line_style);

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
