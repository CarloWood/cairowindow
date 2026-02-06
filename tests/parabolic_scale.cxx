#include "sys.h"
#include "math/QuadraticPolynomial.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/Direction.h"
#include "cairowindow/symbolic/symbolic.h"
#include "utils/almost_equal.h"
#include "utils/square.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include <Eigen/Dense>
#include <sstream>
#include <thread>
#include <ranges>
#include <vector>
#include <array>
#include <random>
#include "debug.h"

//============================================================================
// Code to brute force generate two parabola and the, up till four, points
// at which they differ 10% of the height of the second parabola.
//

constexpr int a = 0;
constexpr int b = 1;
constexpr int c = 2;
constexpr int vx = 3;
constexpr int vy = 4;

double random_double(double min, double max, std::mt19937& gen)
{
  std::uniform_real_distribution<> dis(min, max);
  return dis(gen);
}

// Function to generate random coefficients for a parabola.
std::array<double, 5> generate_parabola(std::mt19937& gen, int given_sign)
{
  std::array<double, 5> P;
  // Generate vertex coordinates.
  P[vx] = random_double(-1.0, 1.0, gen);
  P[vy] = random_double(-1.0, 1.0, gen);
  // Generate c.
  P[c] = random_double(0.02, 0.05, gen);
  if (given_sign == -1 || (given_sign == 0 && random_double(0.0, 1.0, gen) < 0.5))
    P[c] *= -1.0;
  // Calculate corresponding values of a and b.
  P[b] = -2 * P[c] * P[vx];
  P[a] = P[vy] - P[b] * P[vx] - P[c] * P[vx] * P[vx];

  return P;
}

double evaluate(std::array<double, 5> const& P, double x)
{
  return P[a] + (P[b] + P[c] * x) * x;
}

double binary_search(std::array<double, 5> const& P1, std::array<double, 5> const& P2, double left, double right)
{
  constexpr int l = 0;
  constexpr int r = 1;

  std::array<double, 2> x{{left, right}};
  std::array<double, 2> y1;
  std::array<double, 2> y2;
  std::array<int, 2> toggle;

  for (int lr = l; lr <= r; ++lr)
  {
    y1[lr] = evaluate(P1, x[lr]);
    y2[lr] = evaluate(P2, x[lr]);
    toggle[lr] = (std::abs(y1[lr] - y2[lr]) <= 0.1 * std::abs(y2[lr] - P2[vy])) ? -1 : 1;
  }

  for (;;)
  {
    double nx = 0.5 * (x[l] + x[r]);
    double y1t = evaluate(P1, nx);
    double y2t = evaluate(P2, nx);
    int tt = (std::abs(y1t - y2t) <= 0.1 * std::abs(y2t - P2[vy])) ? -1 : 1;

    int lr = (tt == toggle[l]) ? l : r;
    x[lr] = nx;
    y1[lr] = y1t;
    y2[lr] = y2t;

    if (std::abs(x[lr] - x[1 - lr]) < 1e-8 * std::abs(x[lr] + x[1 - lr]))
    {
//      double diff = std::abs(y1t - y2t);
//      double threshold = 0.1 * std::abs(y2t - P2[vy]);
      return nx;
//      std::cout << "x = " << nx << "; y1 = " << y1t << "; y2 = " << y2t << "; diff = " << diff <<
//        "; h2 = " << (y2t - P2[vy]) <<  "; threshold = " << threshold << std::endl;
//      break;
    }
  }
}

// Function to find x coordinates where the absolute difference between two parabolas equals one-tenth of the height of the second parabola.
std::vector<double> find_x_coordinates(std::array<double, 5> const& P1, std::array<double, 5> const& P2)
{
  std::vector<double> xs;
  // Iterate over a range of x coordinates to find the points where the absolute difference meets the condition.
  int toggle = 0;
  for (double x = -100.0; x <= 100.0; x += 0.001)
  {
    double y1 = evaluate(P1, x);
    double y2 = evaluate(P2, x);
    double diff = std::abs(y1 - y2);
    double threshold = 0.1 * std::abs(y2 - P2[vy]);
    if (toggle == 0)
      toggle = (diff <= threshold) ? -1 : 1;
    else if (toggle != ((diff <= threshold) ? -1 : 1))
    {
      xs.push_back(binary_search(P1, P2, x - 0.001, x));
      toggle = -toggle;
    }
  }
  return xs;
}

void print_parabola(std::array<double, 5> const& P, int n)
{
  std::cout << "Parabola " << n << " coefficients: ";
  for (int i = 0; i < 3; ++i)
    std::cout << P[i] << " ";
  std::cout << "\nVertex: " << P[3] << ", " << P[4] << std::endl;
}

std::vector<double> generate_brute_force(std::mt19937& gen, std::array<math::QuadraticPolynomial<double>, 2>& out_parabolas)
{
  for (;;)
  {
    // Generate two parabolas.
    auto P1 = generate_parabola(gen, 0);
    auto P2 = generate_parabola(gen, P1[c] < 0 ? -1 : 1);

    auto xs = find_x_coordinates(P1, P2);

    bool ok = !xs.empty();
    for (double x : xs)
    {
      if (x < -10.0 || x > 10.0)
      {
        ok = false;
        break;
      }
    }
    if (ok)
    {
      out_parabolas[0][0] = P1[a];
      out_parabolas[0][1] = P1[b];
      out_parabolas[0][2] = P1[c];
      out_parabolas[1][0] = P2[a];
      out_parabolas[1][1] = P2[b];
      out_parabolas[1][2] = P2[c];

#ifdef CWDEBUG
      // Display the generated coefficients.
      for (int i = 0; i < 2; ++i)
        std::cout << out_parabolas[i] << std::endl;
#endif

      return xs;
    }
  }
}

// End brute force code.
//============================================================================

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  unsigned seed = 123;
  static std::random_device rd;
  std::mt19937 gen(rd()); //(seed);

  double const w_min = -10.0;
  double const w_max = 10.0;
  int const steps = 100;

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Parabolic scale test", 1200, 900);

    // Create a new layer with a white background.
    auto background_layer = window.create_background_layer<Layer>(color::lightgray COMMA_DEBUG_ONLY("background_layer"));

    // Create another layer.
    auto second_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));

    // Open the window and start drawing.
    std::thread event_loop([&](){
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    });

    double L_min = -6.0;
    double L_max = 6.0;

    // Create and draw plot area.
    plot::Plot plot(window.geometry(), draw::PlotAreaStyle({.color = color::orange}),
        "Parabolic scale test", {},
        "w", {},
        "L", {});
    plot.set_xrange({w_min, w_max});
    plot.set_yrange({L_min, L_max});
    plot.add_to(background_layer, false);

    draw::PointStyle point_style({.color_index = 31, .filled_shape = 1});
    draw::LineStyle line_style({.line_width = 1.0});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::LineStyle derivative_line_style({.line_color = color::turquoise, .line_width = 1.0});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    std::array<draw::LineStyle, 2> parabola_line_styles = {
      draw::LineStyle{{.line_color = color::red, .line_width = 1.0}},
      draw::LineStyle{{.line_color = color::blue, .line_width = 1.0}}
    };

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      uint_fast32_t seed = std::uniform_int_distribution<>{0, 1000000}(gen);
//      seed = 453236;
//      seed = 162473;
//      seed = 283654;
//      seed = 590530;
//      seed = 665162;
//      seed = 969516;
      Dout(dc::notice, "seed = " << seed);
      std::mt19937 gen2(seed);

      // Generate two random parabola's.
      std::array<math::QuadraticPolynomial<double>, 2> parabolas;
      std::vector<double> xs = generate_brute_force(gen2, parabolas);

      // parabolas[0] is the "old" parabola, and parabolas[1] is the "new" parabola.
      math::QuadraticPolynomial<double> const parabola_(parabolas[0]);
      math::QuadraticPolynomial<double> const parabola(parabolas[1]);
      double const v_y = parabola.vertex_y();

      // Draw both parabola, but offset them relative the y-coordinate of the second parabola.
      std::array<plot::BezierFitter, 3> plot_fitters;
      for (int i = 0; i < 2; ++i)
      {
        plot_fitters[i].solve([&](double w) -> Point { return {w, parabolas[i](w) - v_y}; }, plot.viewport());
        plot.add_bezier_fitter(second_layer, parabola_line_styles[i], plot_fitters[i]);
      }
      // If the second parabola is negative, also draw the absolute value.
      if (parabola[2] < 0.0)
      {
        plot_fitters[2].solve([&](double w) -> Point { return {w, v_y - parabolas[1](w)}; }, plot.viewport());
        plot.add_bezier_fitter(second_layer, parabola_line_styles[1], plot_fitters[2]);
      }

      // Draw a vertical line where the distance between the two parabola equals 10% of the height of the second parabola.
      std::array<plot::Line, 4> plot_vertical_lines;
      for (int i = 0; i < xs.size(); ++i)
      {
        plot_vertical_lines[i] = plot::Line{{xs[i], 0.0}, Direction::up};
        plot.add_line(second_layer, line_style, plot_vertical_lines[i]);
      }

      // Draw the difference between the two parabola.
      BezierFitter diff_fitter;
      diff_fitter.solve([&](double w) -> Point{
          auto diff = parabola - parabola_;
          return {w, std::abs(10.0 * diff(w))};
        }, plot.viewport());
      auto plot_diff = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::yellow}), std::move(diff_fitter));

      //----------------------------------------------------------------------
      // Start of algorithm under tests.

      std::array<plot::Line, 4> plot_vertical_lines_algo;

      int count;
      std::array<double, 4> toggles;
      bool close_at_inf = parabola.equal_intervals(parabola_, toggles, count);
      Dout(dc::notice, "close_at_inf = " << std::boolalpha << close_at_inf);
      for (int i = 0; i < count; ++i)
      {
        Dout(dc::notice, "  toggles[" << i << "] = " << toggles[i]);
        plot_vertical_lines_algo[i] = plot::Line{{toggles[i], 0.0}, Direction::up};
        plot.add_line(second_layer,
            line_style({.line_color = color::green, .line_width = 2.0, .dashes = {10.0, 10.0}}), plot_vertical_lines_algo[i]);
      }

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

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
