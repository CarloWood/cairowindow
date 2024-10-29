#include "sys.h"
#include "math/CubicPolynomial.h"
#include "math/AnalyzedCubic.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/almost_equal.h"
#include <random>
#include <chrono>
#include <sstream>
#include <array>
#include "debug.h"

constexpr int number_of_cubics = 1000000;

bool stop = false;

// THIS IS CODE DUPLICATION FROM CubicPolynomial::get_roots!
int get_roots(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer>& layer,
    cairowindow::draw::LineStyle& line_style, cairowindow::plot::BezierFitter& plot_fitter,
    math::CubicPolynomial& cubic, std::array<double, 3>& roots_out, int& iterations)
{
  DoutEntering(dc::notice, "get_roots() for " << cubic);

  double& coefficients_0 = cubic[0];
  std::array<double, 4>& coefficients_ = reinterpret_cast<std::array<double, 4>&>(coefficients_0);

  if (coefficients_[3] == 0.0)
  {
    // The cubic is actually a quadratic.
    math::QuadraticPolynomial qp(coefficients_[0], coefficients_[1], coefficients_[2]);
    return qp.get_roots(*reinterpret_cast<std::array<double, 2>*>(&roots_out[0]));
  }

  // Step one: divide all coefficients by coefficients_[3]. This does not change the roots.
  double const c0 = coefficients_[0] / coefficients_[3];
  double const c1 = coefficients_[1] / coefficients_[3];
  double const c2 = coefficients_[2] / coefficients_[3];
  double const d = utils::square(c2) - 3.0 * c1;

  // The cubic is now monic:
  //
  //   p(x) = c0 + c1 x + c2 x² + x³
  //
  // The first derivative is,
  //
  //   p'(x) = c1 + 2 c2 x + 3 x²
  //
  // Setting this to zero gives:
  //                                  -c2 +/- sqrt(c2^2 - 3c1)   -c2 +/- sqrt(d)
  //   0 = c1 + 2 c2 x + 3 x² --> x = ------------------------ = ---------------
  //                                             3                     3
  // Remember if we have local extrema or not.
  bool const cubic_has_local_extrema = d > 0.0;

  // Calculate the inflection point (where the second derivative is zero).
  double const Ix = -c2 / 3.0;
  // Magic (took me several days to find this).
  double const M = 27.0 * c0 + (2.0 * d - 3.0 * c1) * c2;

  double C0, C1;
  double scale;

  if (cubic_has_local_extrema)
  {
    // Transform the cubic into
    //
    //   Q(u) = C0 - 3u + u³
    //
    // After transforming the cubic we have the following four possibilities:
    //
    //          starting point
    //                ↓
    //  --------------+--O--> u-axis (for example)
    //      .-.         /
    //     /   \       /
    //    /     \←----/----- Inflection point
    //   /       \   /
    //            `-´←------ The extreme with the largest absolute y value.
    //       ↑  ↑  ↑  ↑
    //      -1  0  1  2
    //
    // The starting point is set on the positive slope on the side of the local extreme that is the furthest away from the x-axis.
    //
    double const sqrt_d = std::sqrt(d);
    C0 = M / (d * sqrt_d);
    C1 = -3.0;

    // The applied transform means that any root found must be scaled back by multiplying with
    scale = sqrt_d / 3.0;
    // and then adding Ix back.
  }
  else
  {
    // Transform the cubic into
    //
    //   Q(u) = C0 + 3u + u³
    //
    //          starting point
    //                ↓
    //                /
    //  -------------O+-> u-axis (for example)
    //              /
    //             /
    //            /
    //         .⋅´←--------- Inflection point
    //        /
    //       /
    //      /
    //       ↑  ↑  ↑  ↑
    //      -1  0  1  2

    double const sqrt_md = std::sqrt(-d);
    C0 = M / (-d * sqrt_md);
    C1 = 3.0;

    // The applied transform means that any root found must be scaled back by multiplying with
    scale = sqrt_md / 3.0;
    // and then adding Ix back.
  }

  // Determine if the inflection point is above or below the x-axis.
  bool const inflection_point_y_larger_than_zero = C0 > 0.0;

  // Special case for if zero is a root (i.e. Q(u) = u⋅(u² ± 3)).
  double u = 0.0;

  // Plot the cubic.
  plot_fitter.solve(
      [=](double u) -> cairowindow::Point { return {u, C0 + u * (C1 + utils::square(u))}; },
      plot.viewport());
  plot.add_bezier_fitter(layer, line_style({.line_color = cairowindow::color::green}), plot_fitter);

  if (AI_LIKELY(C0 != 0.0))       // Is 0 not a root of the cubic?
  {
    // Avoid the local extrema and the inflection point because the derivative might be zero there too.
    u = inflection_point_y_larger_than_zero ? -2.0 : 2.0;

    int limit = 100;
    double prev_u;
    double step = std::numeric_limits<double>::infinity();
    double prev_step;
    do
    {
      prev_u = u;
      prev_step = step;
      // Calculate Q(u) = C0 + C1 * u + u^3.
      double Q_u = C0 + u * (utils::square(u) + C1);
      // Calculate Q''(u) = 6 * u;
      double half_Qpp_u = 3.0 * u;
      // Calculate Q'(u) = C1 + 3 * u^2.
      double Qp_u = half_Qpp_u * u + C1;
      // Apply Halley's method.
      step = -Q_u * Qp_u / (utils::square(Qp_u) - Q_u * half_Qpp_u);
      u += step;                                                                // uₙ₊₁ = uₙ - Q(u)Q'(u) / (Q'(u)² - ½Q(u)Q"(u))
//      Dout(dc::notice, "Halley: u = " << std::setprecision(15) << u << "; Δu = " << step);
    }
    while (step != 0.0 /*&& std::abs(step) < std::abs(prev_step)*/ && --limit);
    iterations = 100 - limit;
  }

  Dout(dc::notice, "Root found: " << u);
  stop = std::abs(u) > 15.0;

  roots_out[0] = u * scale + Ix;
  int number_of_roots = 1;

  if (cubic_has_local_extrema)
  {
    // Find the other two roots, if any.
    [[maybe_unused]] double remainder;
    math::QuadraticPolynomial qp = cubic.long_division(roots_out[0], remainder);
    number_of_roots += qp.get_roots(*reinterpret_cast<std::array<double, 2>*>(&roots_out[1]));
  }

  return number_of_roots;
}

void transform(std::array<double, 4>&c, double x_scale, double y_scale, double x_shift, double y_shift)
{
  // First apply the scaling:
  for (int i = 0; i < c.size(); ++i)
    c[i] *= y_scale * std::pow(x_scale, -i);
  // Then the offset.
  c[0] += y_shift + (-c[1] + (c[2] - c[3] * x_shift) * x_shift) * x_shift;
  c[1] += (-2.0 * c[2] + 3.0 * c[3] * x_shift) * x_shift;
  c[2] += -3.0 * c[3] * x_shift;
}

int main(int argc, char* argv[])
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  // Handle random seed.
  unsigned int seed;
  if (argc > 1)
  {
    std::istringstream ss(argv[1]);
    ss >> seed;
  }
  else
    seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::cout << "Seed: " << seed << '\n';

  // Use the mersenne twister engine.
  std::mt19937 engine(seed);

  // Keep all numbers within a reasonable range.
  double const max_value = 1e50;
  double const min_value = 1e-50;
  Dout(dc::notice, "min_value = " << min_value << "; max_value = " << max_value);

  // We'll deal with all numbers logarithmicly, because of their huge spread.
  double const log_max_value = std::log(max_value);
  double const log_min_value = std::log(min_value);
  Dout(dc::notice, "log(min_value) = " << log_min_value << "; log(max_value) = " << log_max_value);

  // Truncate these exponents towards zero, with a tiny little bit more leeway.
  double const log_min_value_trunc = std::trunc(log_min_value) + 1;
  double const log_max_value_trunc = std::trunc(log_max_value) - 1;
  Dout(dc::notice, "log(min_value) = " << log_min_value_trunc << "; log(max_value) = " << log_max_value_trunc);

  std::normal_distribution normal_dist(0.0, log_max_value_trunc / 5.0); // 1 in 3000000 will have the most extreme value.
  auto random_value = [&](){
    double exponent;
    do { exponent = normal_dist(engine); } while (exponent < log_min_value_trunc || log_max_value_trunc > log_max_value_trunc);
    return std::exp(exponent);
  };

  // Generate random cubics.
  std::cout << "Generating random cubic polynomials..." << std::endl;
  std::vector<math::CubicPolynomial> cubics;
  for (int i = 0; i < number_of_cubics; ++i)
  {
    // Generate the position of the inflection point.
    double Ix = random_value();
    double Iy = random_value();
    // Generate the position of the local minimum.
    double Ex = random_value();
    double Ey = random_value();
    // Create an array with coefficients for the cubic x * (x^2 - 3) = -3 x + x^3 if Ey < Iy, or +3 x + x^3 if Ey >= Iy.
    std::array<double, 4> c = { 0, std::copysign(3.0, Ey - Iy), 0, 1 };
    // Scale the cubic such that the minimum of -3 x + x^3 (at x = 1) goes to Ex - Ix:
    double x_scale = Ex - Ix;
    // Scale the cubic such that the minimum, at y = -2, goes to Ey - Iy.
    double y_scale = 0.5 * (Iy - Ey);
    // Shift the cubic such that the inflection point (still at the origin) moves to (Ix, Iy).
    double x_shift = Ix;
    double y_shift = Iy;

    // Apply the linear transformation.
    transform(c, x_scale, y_scale, x_shift, y_shift);

    cubics.emplace_back(c[0], c[1], c[2], c[3]);
  }

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("random_cubics", 1600, 1200);

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
        "Transformed cubic to calculate roots", {},
        "u", {},
        "v", {});
    plot.set_xrange({-15.0, 15.0});
    plot.set_yrange({-100.0, 100.0});
    plot.add_to(background_layer, false);

    draw::LineStyle curve_line_style({.line_width = 1.0});

    std::array<double, 3> roots;
    std::cout << "Calculating roots..." << std::endl;
    int max_iterations = 0;
//    Debug(libcw_do.off());
    unsigned long total_iterations = 0;
    for (math::CubicPolynomial& cubic : cubics)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      int iterations;
      plot::BezierFitter plot_bezier_fitter;
      int n = get_roots(plot, second_layer, curve_line_style, plot_bezier_fitter, cubic, roots, iterations);
      total_iterations += iterations;

      if (iterations == 12 || (iterations < 100 && iterations > max_iterations))
      {
        max_iterations = iterations;
//        Debug(libcw_do.on());

        bool inflection_point_y_larger_than_zero =
          13.5 * cubic[0] + cubic[2] * (utils::square(cubic[2] / cubic[3]) - 4.5 * (cubic[1] / cubic[3])) > 0.0;

        math::AnalyzedCubic acubic;
        // Calculate the minimum (-1) if the inflection point lays above the x-axis, and the maximum otherwise.
        // This way acubic.get_extreme() (E above) becomes the extreme that is the closest to the x-axis.
        acubic.initialize(cubic, inflection_point_y_larger_than_zero ? -1 : 1);

        // Obtain the calculated inflection point.
        double const inflection_point_x = acubic.inflection_point();

        // Remember if we have local extrema or not.
        bool const cubic_has_local_extrema = acubic.has_extrema();

        // Avoid the local extrema and the inflection point because the derivative might be zero there too.
        double x = cubic_has_local_extrema ?
            3 * inflection_point_x - 2 * acubic.get_extreme() :
            inflection_point_x + (((cubic[3] > 0.0) == inflection_point_y_larger_than_zero) ? -1.0 : 1.0);

        Dout(dc::notice|continued_cf, cubic << " roots: ");
        for (int r = 0; r < n; ++r)
          Dout(dc::continued, roots[r] << ", ");
        Dout(dc::continued, "inflection point: " << inflection_point_x);
        if (cubic_has_local_extrema)
          Dout(dc::continued, ", extreme: " << acubic.get_extreme());
        Dout(dc::finish, ", starting point: " << x);
//        Debug(libcw_do.off());
      }

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      if (stop)
      {
        // Block until a redraw is necessary (for example because the user moved a draggable object,
        // or wants to print the current drawing) then go to the top of loop for a redraw.
        if (!window.handle_input_events())
          break;          // Program must be terminated.
      }
    }
//    Debug(libcw_do.on());
    Dout(dc::notice, "Average number of iterations per cubic: " << (total_iterations / number_of_cubics));

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
