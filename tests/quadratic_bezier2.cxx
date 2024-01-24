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

#define USE_P_BETA 0
#define USE_P_GAMMA 0

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
    draw::TextStyle<> slider_style{.position = draw::centered_below, .font_size = 18.0, .offset = 10};
    draw::LineStyle curve_line_style{.line_width = 1.0};
    draw::LineStyle solid_line_style{.line_color = color::black, .line_width = 1.0};
    draw::LineStyle line_style{.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}};
    draw::ArcStyle arc_style{.line_color = color::blue, .line_width = 1.0};

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

    // Create a point Q on the circle.
    double phi = M_PI;
    auto plot_Q = plot.create_point(second_layer, P0P1_circle_center + 0.5 * Direction{phi}, point_style);

#if USE_P_BETA
    // Initial position of Pᵦ, a point at t < 0.
    auto plot_P_beta = plot.create_point(second_layer, {-0.5, -0.25}, point_style);
#if USE_P_GAMMA
    // Initial position of Pᵧ, a point at t > 1.
    auto plot_P_gamma = plot.create_point(second_layer, {1.8, -0.15}, point_style);
#endif
#endif

    // Draw a slider for w.
    auto slider_w = plot.create_slider(second_layer, {928, 83, 7, 400}, 0.5, 0.25, 2.0);
    auto slider_w_label = plot.create_text(second_layer, Pixel{928, 483}, "w", slider_style);

    // Draw a slider for s.
    auto slider_s = plot.create_slider(second_layer, {978, 83, 7, 400}, 2.0, 0.5, 4.0);
    auto slider_s_label = plot.create_text(second_layer, Pixel{978, 483}, "s", slider_style);

    // Draw a slider for v.
    auto slider_v = plot.create_slider(second_layer, {1028, 83, 7, 400}, 0.5, -2.0, 2.0);
    auto slider_v_label = plot.create_text(second_layer, Pixel{1028, 483}, "v", slider_style);

    // Allow dragging Pᵦ and Pᵧ.
#if USE_P_BETA
    window.register_draggable(plot, &plot_P_beta);
#if USE_P_GAMMA
    window.register_draggable(plot, &plot_P_gamma);
#endif
#endif
    window.register_draggable(plot, &plot_Q, [&P0P1_circle_center](Point const& new_position)
        {
          return P0P1_circle_center + 0.5 * Direction{P0P1_circle_center, new_position};
        }
    );

    double prev_w = 0.0;
    double prev_s{}, prev_v{};
    Point prev_Q = plot_Q;
    bool w_is_negative = false;

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

#if USE_P_BETA
      // Draw a label for Pᵦ.
      auto P_beta_label = plot.create_text(second_layer, plot_P_beta, "Pᵦ", label_style({.position = draw::centered_right_of}));

#if USE_P_GAMMA
      // Draw a label for Pᵧ.
      auto P_gamma_label = plot.create_text(second_layer, plot_P_gamma, "Pᵧ", label_style({.position = draw::centered_right_of}));
#endif
#endif

      // Draw a label for Q.
      auto Q_label = plot.create_text(second_layer, plot_Q, "Q", label_style);

      // Draw a line from P₀ to Q.
      auto plot_P0_Q = plot.create_connector(second_layer, plot_P0, plot_Q, line_style);
      Vector P0_Q{plot_P0, plot_Q};
      auto P0_Q_label = plot.create_text(second_layer, plot_P0 + 0.5 * P0_Q, "s²(1-2v)",
          label_style({.font_size = 12, .offset = 5}));

      // Draw a line between P₁ and Q.
      auto line_P1_Q = plot.create_connector(second_layer,
          P1, plot_Q, Connector::open_arrow, Connector::open_arrow, line_style({.line_color = color::coral, .dashes = {3.0, 3.0}}));
      Vector Q_P1{plot_Q, P1};
      auto sw_label2 = plot.create_text(second_layer, plot_Q + 0.5 * Q_P1,
          "sw", label_style({.position = draw::centered_below, .font_size = 12, .offset = 5}));

      // Determine the distance between Q and respectively P₀ and P₁.
      double s_times_w = Q_P1.length();
      Direction perpendicular_to_symmetry_line_dir = (s_times_w > 1e-6) ? Q_P1.direction() : Direction::down;
      if (w_is_negative)
        perpendicular_to_symmetry_line_dir = perpendicular_to_symmetry_line_dir.inverse();
      // Calculate s_squared_times_one_minus_two_v assuming w_is_negative won't change.
      Direction symmetry_line_dir = perpendicular_to_symmetry_line_dir.normal();
      double s_squared_times_one_minus_two_v = P0_Q.dot(symmetry_line_dir);

      Dout(dc::notice, "s²(1-2v) = " << s_squared_times_one_minus_two_v);
      Dout(dc::notice, "sw = " << (w_is_negative ? -s_times_w : s_times_w));

      // Read-out the slider values;
      double w = slider_w.value();
      double s = slider_s.value();
      double v = slider_v.value();
      double Qx = plot_Q.x();
      double Qy = plot_Q.y();

      static constexpr double min_s = 0.1;
      bool slider_values_changed = false;
      if (w != prev_w || prev_Q != plot_Q)
      {
        if (((plot_Q.y() >= 0.0) != (prev_Q.y() >= 0.0)) && plot_Q.x() > 0.5)
        {
          w_is_negative = !w_is_negative;
          perpendicular_to_symmetry_line_dir = perpendicular_to_symmetry_line_dir.inverse();
          symmetry_line_dir = symmetry_line_dir.inverse();
          s_squared_times_one_minus_two_v = -s_squared_times_one_minus_two_v;
        }
        prev_Q = plot_Q;

        // Calculate the values of s and v.
        s = std::max(min_s, s_times_w / w);    // s must always be larger than zero.
        v = (1.0 - s_squared_times_one_minus_two_v / (s * s)) / 2.0;
        slider_values_changed = true;
      }
      else if (s != prev_s)
      {
        // Calculate the values of w and v.
        w = s_times_w / s;
        v = (1.0 - s_squared_times_one_minus_two_v / (s * s)) / 2.0;
        slider_values_changed = true;
      }
      else if (v != prev_v)
      {
        // Calculate the values of s and w.
        double s_squared = s_squared_times_one_minus_two_v / (1.0 - 2.0 * v);
        Dout(dc::notice, "v = " << v << "; s_squared = " << s_squared);
        if (s_squared < min_s * min_s)
        {
          // This means we need to limit v to some value near 0.5.
          // If s_squared became negative then that is because v passed the value of 0.5,
          // in other words, s went to +inf - and we should limit it to slider_s.max_value().
          // However if s_squared is still positive then it is just too small and we want
          // to limit it to slider_s.min_value().
          s_squared = (s_squared < 0) ? slider_s.max_value() * slider_s.max_value() : slider_s.min_value() * slider_s.min_value();
          v = (1.0 - s_squared_times_one_minus_two_v / s_squared) / 2.0;
          Dout(dc::notice, "Changed v to " << v);
        }
        s = std::sqrt(s_squared_times_one_minus_two_v / (1.0 - 2.0 * v));
        w = s_times_w / s;
        Dout(dc::notice, "s = " << s << "; w = " << w);
        slider_values_changed = true;
      }

      // If any slider changed, set all sliders to the new (calculated) values,
      // and read them out again in order to be sure that we can use != on the
      // double values to detect if they were changed, the next time we reenter
      // this loop.
      if (slider_values_changed)
      {
        slider_w.set_value(std::max(0.0001, w));        // w must always be larger than zero.
        slider_s.set_value(s);
        slider_v.set_value(v);

        prev_w = slider_w.value();
        prev_s = slider_s.value();
        prev_v = slider_v.value();
      }

      Dout(dc::notice, "w = " << w << "; s = " << s << "; v = " << v);

#if USE_P_BETA
      double x_beta = plot_P_beta.x();
      double y_beta = plot_P_beta.y();
#if USE_P_GAMMA
      double x_gamma = plot_P_gamma.x();
      double y_gamma = plot_P_gamma.y();

      // Helper variables.
      double x_span = x_gamma - x_beta;
      double y_span = y_gamma - y_beta;
      double subexpr50 = y_span * ((1.0 - x_gamma) * x_gamma / y_gamma - (1.0 - x_beta) * x_beta / y_beta);
      double tan_theta = (std::sqrt(x_span * x_span + subexpr50) - x_span) / y_span;
      double z = x_beta + y_beta * tan_theta;

      // Define the matrix M.
      double m10 = y_beta / (z * (1.0 - z));
      double m01 = m10 * tan_theta;
      double m11 = -m10;
      double m00 = 1.0 - m01;
#else
#endif
#else
      // Helper variables.
      double theta = perpendicular_to_symmetry_line_dir.as_angle();
      double sw = s * w;
      if (w_is_negative)
        sw = -sw;

      // Define the matrix M.
      double m00 = sw * std::cos(theta) + 2.0 * s * s * v * std::sin(theta);
      double m01 = -s * s * std::sin(theta);
      double m10 = sw * std::sin(theta) - 2.0 * s * s * v * std::cos(theta);
      double m11 = s * s * std::cos(theta);
#endif

      Dout(dc::notice, "M = ((" << m00 << ", " << m01 << "), (" << m10 << ", " << m11 << "))");

      if (std::isnan(m10))
      {
        // Can't draw the curve.
        window.set_send_expose_events(true);
        window.handle_dragging();
        continue;
      }

      auto xt = [=](double t){ return t * (m00 + m01 * t); };
      auto yt = [=](double t){ return t * (m10 + m11 * t); };

      std::vector<Point> curve_points;
      for (int i = -400; i <= 400; ++i)
      {
        double t = i * 0.01;
        Vector v{xt(t), yt(t)};
        curve_points.push_back(v.point());
      }
      auto curve = plot.create_curve(second_layer, std::move(curve_points), curve_line_style);

      // V, the parabola vertex point resides at t=v.
      auto V = plot.create_point(second_layer, {xt(v), yt(v)}, point_style);
      label_style.position = draw::centered_left_of;
      auto V_label = plot.create_text(second_layer, V, "V", label_style({.position = draw::centered_below}));

      // Draw a vertical line from V up 0.5.
      Point V6(V.x(), V.y() + 0.5);
      LinePiece foo(V, V6);
      auto plot_foo = plot.create_line(second_layer, foo, line_style);

      // Point on symmetry line at distance 1 from V.
      Point V1 = V + symmetry_line_dir;
      auto plot_V1 = plot.create_point(second_layer, V1, point_style);
      // Draw an arrow from V to V1.
      auto plot_V_V1 = plot.create_connector(second_layer, V, V1, solid_line_style({.line_color = color::purple}));

      // Draw the unit vectors X and Y.
      auto plot_X = plot.create_connector(second_layer, V, V + perpendicular_to_symmetry_line_dir, solid_line_style);
      auto X_label = plot.create_text(second_layer, V + 0.5 * perpendicular_to_symmetry_line_dir, "X", label_style);
      auto plot_Y = plot.create_connector(second_layer, V, V + symmetry_line_dir, solid_line_style);
      auto Y_label = plot.create_text(second_layer, V + 0.5 * symmetry_line_dir, "Y", label_style);

      // Go from V1 a distance w left and right.
      Point V1L = V1 + w * perpendicular_to_symmetry_line_dir;
      Point V1R = V1 - w * perpendicular_to_symmetry_line_dir;
      auto plot_V1L = plot.create_point(second_layer, V1L, point_style({.color_index = color_index2}));
      auto plot_V1R = plot.create_point(second_layer, V1R, point_style({.color_index = color_index2}));

      // Draw the symmetry line of the parabola, through V.
      auto symmetry_line_of_parabola = plot.create_line(second_layer, V, symmetry_line_dir, line_style);

      // Draw a line between V1 and V1L.
      auto line_at_one_from_V = plot.create_connector(second_layer,
          V1, V1L, Connector::open_arrow, Connector::open_arrow, line_style({.line_color = color::coral, .dashes = {3.0, 3.0}}));
      auto w_label = plot.create_text(second_layer, V1 + 0.5 * w * perpendicular_to_symmetry_line_dir,
          "w", label_style({.position = draw::centered_above, .font_size = 14, .offset = 5}));

      // Draw an arc between P1_Q and P1_P0.
      plot::Arc alpha_arc(P1, perpendicular_to_symmetry_line_dir.inverse(), Direction{P1, plot_P0}, 0.1);
      plot.add_arc(second_layer, alpha_arc, arc_style({.line_color = color::blue}));
      auto alpha_label = plot.create_text(second_layer, P1 + 0.13 * alpha_arc.bisector_direction(),
          "α", label_style({.position = draw::centered, .font_size = 14}));

      plot::Arc theta_arc(V, Direction::up, symmetry_line_dir, 0.2);
      plot.add_arc(second_layer, theta_arc, arc_style);
      auto theta_label = plot.create_text(second_layer, V + 0.24 * theta_arc.bisector_direction(),
          "θ", label_style({.position = draw::centered, .font_size = 14}));

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
