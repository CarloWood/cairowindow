#include "sys.h"
#include "CubicPolynomial.h"
#include "gradient_descent/Sample.h"
#include "gradient_descent/HorizontalDirection.h"
#include "gradient_descent/VerticalDirection.h"
#include "gradient_descent/Approximation.h"
#include "gradient_descent/Weight.h"
#include "cairowindow/BezierCurve.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/draw/Point.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/ColorPool.h"
#include "utils/square.h"
#include <random>
#include "debug.h"

bool check_hdirection(gradient_descent::HorizontalDirection hdirection, double w, double wl, double wr)
{
  ASSERT(wl < wr);
  if (w < wl)
    return hdirection == gradient_descent::HorizontalDirection::left;
  else if (w > wr)
    return hdirection == gradient_descent::HorizontalDirection::right;
  return hdirection == gradient_descent::HorizontalDirection::inbetween;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

// std::random_device rd;
  int seed = 10248;
  Dout(dc::notice, "seed = " << seed);
  std::mt19937 engine(seed);

  std::uniform_real_distribution<double> a_dist(-200.0, 200.0);
  std::uniform_real_distribution<double> b_dist(-200.0, 200.0);
  std::uniform_real_distribution<double> c_dist(-20.0, 20.0);
  std::uniform_real_distribution<double> d_dist(-1.0, 1.0);

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;
    using HorizontalDirection = gradient_descent::HorizontalDirection;
    using VerticalDirection = gradient_descent::VerticalDirection;

    // Create a window.
    Window window("Cubic", 1200, 900);

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

    utils::ColorPool<32> color_pool;
    draw::PointStyle point_style({.color_index = color_pool.get_and_use_color(), .filled_shape = 1});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::LineStyle solid_line_style({.line_color = color::black, .line_width = 1.0});
    draw::LineStyle line_style({.line_color = color::black, .line_width = 1.0, .dashes = {10.0, 5.0}});
    draw::ArcStyle arc_style({.line_color = color::blue, .line_width = 1.0});

    double wl;
    double wr;
    double ll;
    double rl;
    math::CubicPolynomial cubic;
    gradient_descent::Sample sl;
    gradient_descent::Sample sr;
    double minimum;
    double maximum;
    double D;

    while (true)
    {
      double a = a_dist(engine);
      double b = b_dist(engine);
      double c = c_dist(engine);
      double d = d_dist(engine);

      cubic[0] = a;
      cubic[1] = b;
      cubic[2] = c;
      cubic[3] = d;

      Dout(dc::notice, "cubic = " << cubic);

      D = utils::square(c) - 3.0 * b * d;

      if (D <= 0.0)
      {
        // Find the point where the first derivative is minimal.
        // f(x) = a + b x + c x^2 + d x^3
        // f'(x) = b + 2c x + 3d x^2
        // f''(x) = 2c + 6d x = 0 --> x = -c/3d
        double x = -c / (3.0 * d);
        // The value of the derivative at that point:
        double dfdw_min = b + (2.0 * c + 3.0 * d * x) * x;
        double dfdw_aim;
        if (std::abs(dfdw_min) < 0.1)
          dfdw_aim = (d < 0.0) ? -1.0 : 1.0;
        else
          dfdw_aim = 10.0 * dfdw_min;
        // Solve,
        // f'(x) = b + 2c x + 3d x^2 = dfdw_aim -->
        // 3d x^2 + 2c x + b - dfdw_aim = 0 -->
        // x = (-2c +/- sqrt(4c^2 - 12d (b - dfdw_aim))) / (6d)
        double D2 = 4.0 * c * c - 12.0 * d * (b - dfdw_aim);
        ASSERT(D2 > 0.0);
        ll = (-2.0 * c - std::sqrt(D2)) / (6.0 * d);
        rl = (-2.0 * c + std::sqrt(D2)) / (6.0 * d);

        std::uniform_real_distribution<double> wrl_dist(ll, rl);
        wl = wrl_dist(engine);
        do
        {
          wr = wrl_dist(engine);
        }
        while (std::abs(wl - wr) < 1e-3);
        if (wl > wr)
          std::swap(wl, wr);
      }
      else
      {
        double sqrt_D = std::sqrt(D);
        minimum = (-c + sqrt_D) / (3.0 * d);
        maximum = (-c - sqrt_D) / (3.0 * d);
        Dout(dc::notice, "minimum = " << minimum << ", maximum = " << maximum);
        //          lm         rm
        // ll-------max--------min--------rl  d>0
        // ll-------min--------max--------rl d<0
        // 0  wl,wr
        // 1  wl          wr
        // 2  wl                     wr
        // 3            wl,wr
        // 4              wl         wr
        // 5                       wl,wr
        std::uniform_int_distribution<int> case_dist(0, 5);
        int c = case_dist(engine);
        double lm = std::min(minimum, maximum);
        double rm = std::max(minimum, maximum);
        ll = lm - (rm - lm);
        rl = rm + (rm - lm);
        std::uniform_real_distribution<double> left_dist(ll, lm);
        std::uniform_real_distribution<double> inbetween_dist(lm, rm);
        std::uniform_real_distribution<double> right_dist(rm, rl);
        switch (c)
        {
          case 0:
          {
            wl = left_dist(engine);
            do
            {
              wr = left_dist(engine);
            }
            while (std::abs(wl - wr) < 1e-3);
            if (wl > wr)
              std::swap(wl, wr);
            break;
          }
          case 1:
          {
            wl = left_dist(engine);
            wr = inbetween_dist(engine);
            break;
          }
          case 2:
          {
            wl = left_dist(engine);
            wr = right_dist(engine);
            break;
          }
          case 3:
          {
            wl = inbetween_dist(engine);
            do
            {
              wr = inbetween_dist(engine);
            }
            while (std::abs(wl - wr) < 1e-3);
            if (wl > wr)
              std::swap(wl, wr);
            break;
          }
          case 4:
          {
            wl = inbetween_dist(engine);
            wr = right_dist(engine);
            break;
          }
          case 5:
          {
            wl = right_dist(engine);
            do
            {
              wr = right_dist(engine);
            }
            while (std::abs(wl - wr) < 1e-3);
            if (wl > wr)
              std::swap(wl, wr);
            break;
          }
        }
        ASSERT(ll < wl && wl < wr && wr < rl);
      }
      sl = gradient_descent::Sample{wl, cubic(wl), cubic.derivative(wl)};
      sr = gradient_descent::Sample{wr, cubic(wr), cubic.derivative(wr)};
      Dout(dc::notice, "sl = " << sl << ", sr = " << sr);

      std::uniform_int_distribution<int> dir_dist(-1, 1);
      HorizontalDirection hdirection = static_cast<HorizontalDirection>(dir_dist(engine));
      VerticalDirection vdirection = static_cast<VerticalDirection>(dir_dist(engine));
      Dout(dc::notice, "hdirection = " << hdirection << ", vdirection = " << vdirection);

      gradient_descent::Approximation approximation;
      approximation.add(&sl, false);
      auto ar = approximation.add(&sr, false);
      if (ar == gradient_descent::ScaleUpdate::first_sample)
        continue;

      HorizontalDirection const input_hdirection = hdirection;
      VerticalDirection const input_vdirection = vdirection;
      gradient_descent::Weight result = approximation.find_extreme(hdirection, vdirection);
      Dout(dc::notice, "result = " << result << "; hdirection = " << hdirection << "; vdirection = " << vdirection);

      // Create and draw plot area.
      plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
          "Clamp Check", {},
          "x", {},
          "y", {});

      {
        double w1 = std::min(wl, ll);
        double w2 = std::max(wr, rl);
        double dw = w2 - w1;
        double const w_min = w1 - 0.2 * dw;
        double const w_max = w2 + 0.2 * dw;

        double Lw1 = std::min(cubic(w1), cubic(w2));
        double Lw2 = std::max(cubic(w1), cubic(w2));
        double dLw = Lw2 - Lw1;
        double const Lw_min = Lw1 - 0.2 * dLw;
        double const Lw_max = Lw2 + 0.2 * dLw;

        plot.set_xrange({w_min, w_max});
        plot.set_yrange({Lw_min, Lw_max});
        plot.add_to(background_layer, false);
      }

      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      // Create a point P₀.
      auto plot_P0 = plot.create_point(second_layer, point_style, {sl.w(), sl.Lw()});
      // Create a point P₁.
      auto plot_P1 = plot.create_point(second_layer, point_style, {sr.w(), sr.Lw()});

      // Draw a label for P₀.
      auto P0_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P0, "P₀");

      // Draw a label for P₁.
      auto P1_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P1, "P₁");

      // Plot the cubic.
      BezierFitter fitter2([&](double w) -> Point{ return {w, cubic(w)}; }, plot.viewport());
      auto plot_cubic = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(fitter2));

      plot::Line plot_minimum_line;
      plot::Line plot_maximum_line;
      if (vdirection != VerticalDirection::unknown)
      {
        ASSERT(D > 0.0);

        if (vdirection == VerticalDirection::down)
        {
          // Draw a vertical line where the minimum is.
          plot_minimum_line = Line{Point{minimum, 0.0}, Direction::up};
          plot.add_line(second_layer, line_style({.line_color = color::red}), plot_minimum_line);
        }
        else if (vdirection == VerticalDirection::up)
        {
          // Draw a vertical line where the maximum is.
          plot_maximum_line = Line{Point{maximum, 0.0}, Direction::up};
          plot.add_line(second_layer, line_style({.line_color = color::red}), plot_maximum_line);
        }
      }
      plot::Line plot_center_line;
      if (input_hdirection == HorizontalDirection::undecided && input_vdirection == VerticalDirection::unknown)
      {
        // Draw a vertical line at the center.
        double center = 0.5 * (wl + wr);
        plot_center_line = Line{Point{center, 0.0}, Direction::up};
        plot.add_line(second_layer, line_style, plot_center_line);
      }

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      Dout(dc::notice, "input_hdirection = " << input_hdirection << ", input_vdirection = " << input_vdirection);

      // Is there no extreme?
      if (D <= 0.0)
      {
        // If there is no extreme, then vdirection should be 'unknown'.
        ASSERT(vdirection == VerticalDirection::unknown);
        // hdirection should be untouched.
        ASSERT(hdirection == input_hdirection);
        continue;
      }

      // Do we not care where the returned extreme is?
      if (input_hdirection == HorizontalDirection::undecided)
      {
        // Then some extreme must be returned.
        ASSERT(vdirection != VerticalDirection::unknown);
        // If we requested a certain type of extreme, then that must be returned.
        ASSERT(input_vdirection == VerticalDirection::unknown || vdirection == input_vdirection);
        // hdirection must be set correctly.
        double w = result;
        ASSERT(check_hdirection(hdirection, w, wl, wr));
        // The result must correspond to the correct extreme.
        double extreme = (vdirection == VerticalDirection::down) ? minimum : maximum;
        ASSERT(std::abs(w - extreme) < 0.1 * std::abs(minimum - maximum));
        continue;
      }

      // Was an extreme was returned?
      if (vdirection != VerticalDirection::unknown)
      {
        // If we requested a certain type of extreme, then that must be returned.
        ASSERT(input_vdirection == VerticalDirection::unknown || vdirection == input_vdirection);
        // hdirection must be set correctly.
        double w = result;
        ASSERT(check_hdirection(hdirection, w, wl, wr));
        // The result must correspond to the correct extreme.
        double extreme = (vdirection == VerticalDirection::down) ? minimum : maximum;
        ASSERT(std::abs(w - extreme) < 0.1 * std::abs(minimum - maximum));
        // The returned hdirection must correspond to the requested one.
        ASSERT(hdirection == HorizontalDirection::inbetween || hdirection == input_hdirection);
        continue;
      }

      // Run over the two possible extremes.
      for (int e = 0; e < 2; ++e)
      {
        double extreme = (e == 0) ? minimum : maximum;
        // Skip extremes that should be rejected because of input_vdirection.
        if ((e == 0 && input_vdirection == VerticalDirection::up) ||
            (e == 1 && input_vdirection == VerticalDirection::down))
          continue;
        // Then this extreme should be rejected because of input_hdirection.
        ASSERT((input_hdirection == HorizontalDirection::left && extreme > std::max(wl, wr)) ||
            (input_hdirection == HorizontalDirection::right && extreme < std::min(wl, wr)));
      }
      // Everything was checked with ASSERTs now.
      continue;

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
