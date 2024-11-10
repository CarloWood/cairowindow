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

#define INTERACTIVE 0

constexpr int number_of_cubics = 1000000;

bool stop = true;

int get_roots(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer>& layer,
    cairowindow::draw::LineStyle& line_style, cairowindow::plot::BezierFitter* plot_fitter,
    math::CubicPolynomial& cubic, std::array<double, 3>& roots_out, double& initial_guess, int& iterations)
{
  DoutEntering(dc::notice, "get_roots() for " << cubic);

  // Define a coefficients_ for use by math/CubicPolynomial_get_roots.cpp.
  std::array<double, 4>& coefficients_ = reinterpret_cast<std::array<double, 4>&>(cubic[0]);

  using math::QuadraticPolynomial;
  // Include the body of the function.
# define RANDOM_CUBICS_TEST
# include "math/CubicPolynomial_get_roots.cpp"
# undef RANDOM_CUBICS_TEST
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

#if 1
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
#endif

#if INTERACTIVE
  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("random_cubics", 1200, 900);
    Window window2("difference", 1200, 900);
    Window window3("S(C0)", 1200, 900);

    // Create a new layer with a gray background.
    auto background_layer = window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));
    auto background_layer2 = window2.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));
    auto background_layer3 = window3.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"));

    // Create another layer.
    auto second_layer = window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));
    auto second_layer2 = window2.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));
    auto second_layer3 = window3.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"));

    // Open the window and start drawing.
    std::thread event_loop([&](){
      Debug(NAMESPACE_DEBUG::init_thread("event_loop"));
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    });

    // Open window2 and start drawing.
    std::thread event_loop2([&](){
      Debug(NAMESPACE_DEBUG::init_thread("event_loop2"));
      EventLoop event_loop = window2.run();
      event_loop.set_cleanly_terminated();
    });

    // Open window3 and start drawing.
    std::thread event_loop3([&](){
      Debug(NAMESPACE_DEBUG::init_thread("event_loop3"));
      EventLoop event_loop = window3.run();
      event_loop.set_cleanly_terminated();
    });

    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "Transformed cubic to calculate roots", {},
        "u", {},
        "C0 - 3u + u^3", {});
    plot.set_xrange({-5.0, 5.0});
    plot.set_yrange({-10.0, 10.0});
    plot.add_to(background_layer, false);

    draw::LineStyle line_style({.line_color = color::red, .line_width = 1.0, .dashes = {10.0, 5.0}});

    auto diff_lambda = [&](double C0) -> cairowindow::Point
    {
      math::CubicPolynomial p(C0, -3, 0, 1);
      std::array<double, 3> roots;
      int iterations;
      double initial_guess;
      int n = get_roots(plot, second_layer, line_style, nullptr, p, roots, initial_guess, iterations);
      return {C0, initial_guess - roots[0]};
    };

    auto rel_err_lambda = [&](double C0) -> cairowindow::Point
    {
      math::CubicPolynomial p(C0, -3, 0, 1);
      std::array<double, 3> roots;
      int iterations;
      double initial_guess;
      int n = get_roots(plot, second_layer, line_style, nullptr, p, roots, initial_guess, iterations);
      return {C0, (initial_guess - roots[0]) / std::abs(roots[0])};
    };

    plot::Plot plot2(window2.geometry(), { .grid = {.color = color::orange} },
        "Relative error of initial guess of root", {},
        "C0", {},
        "(initial_guess - root) / |root|", {});
    plot2.set_xrange({-10.0, 10.0});

    Debug(libcw_do.off());
    double ymin = 1e100;
    double ymax = -1e100;
    for (int i = 1; i <= 100; i += 2)
    {
      double C0 = plot2.xrange().min() + i * plot2.xrange().size() / 100;
      double y = diff_lambda(C0).y();
      ymax = std::max(ymax, y);
      ymin = std::min(ymin, y);
    }
    Debug(libcw_do.on());

    plot2.set_yrange({ymin, ymax});
    plot2.add_to(background_layer2, false);

    plot::Plot plot3(window3.geometry(), { .grid = {.color = color::orange} },
        "Smoothing function", {},
        "C0", {},
        "S(C0)", {});
    plot3.set_xrange({-1.0, 1.0});
    plot3.set_yrange({0.0, 1.0});
    plot3.add_to(background_layer3, false);

    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});

    auto slider_c0 = plot.create_slider(second_layer, {928, 83, 7, 400}, 0.0, plot2.xrange().min(), plot2.xrange().max());
    auto slider_c0_label = plot.create_text(second_layer, slider_style, Pixel{928, 483}, "C0");

    auto slider_s = plot3.create_slider(second_layer3, {928, 83, 7, 400}, 0.0, -1.0, 1.0);
    auto slider_s_label = plot3.create_text(second_layer3, slider_style, Pixel{928, 483}, "s");

    auto SC0_lambda = [&](double C0) -> cairowindow::Point
    {
      math::CubicPolynomial p(C0, -3, 0, 1);
      std::array<double, 3> roots;
      int iterations;
      double initial_guess;
      int n = get_roots(plot, second_layer, line_style, nullptr, p, roots, initial_guess, iterations);
      double root = roots[0];

      //---------------------------------------------------------------
      // Copied from math/CubicPolynomial_get_roots.cpp
      static constexpr double sqrt3 = 1.7320508075688773;    // -√3
      // Define the value of the root in the case that C0 is zero.
      double root0 = std::copysign(sqrt3, -C0);

      // We need the cube root of C0.
      double const cbrtC0 = std::cbrt(C0);

      // Define the value of the root in the case that C0 is large.
      double const root1 = -(cbrtC0 + 1.0 / cbrtC0);
      //---------------------------------------------------------------

      // root = (1 - S(C0)) * root0 + S(C0) * root1 -->
      double SC0 = (root - root0) / (root1 - root0);
      return {C0, SC0};
    };
#endif

#if 0
    //---------------------
    double x_min = -1.0;
    double x_max = 1.0;
    double y_min = SC0_lambda(x_min).y();
    double y_max = SC0_lambda(x_max).y();

    // Check if the solution is bracketed
    if (y_min > 0.5 || y_max < 0.5) {
        throw std::runtime_error("y = 0.5 is not bracketed within x_min and x_max.");
    }

    double tol = 1e-12;
    while ((x_max - x_min) > tol) {
        double x_mid = 0.5 * (x_min + x_max);
        double y_mid = SC0_lambda(x_mid).y();

        if (y_mid < 0.5) {
            x_min = x_mid;
        } else {
            x_max = x_mid;
        }
    }

    Dout(dc::notice, std::setprecision(12) << "Center = " << (0.5 * (x_min + x_max)) << "; aka " << std::pow(10.0, 0.5 * (x_min + x_max)));
    return 0;
    //---------------------
    // Center = -0.0439481247927
#endif
    constexpr double center = -0.0439481247927;

    std::array<double, 3> roots;
    std::cout << "Calculating roots..." << std::endl;
    int max_iterations = 0;
    unsigned long total_iterations = 0;
    int count = 0;
#if INTERACTIVE
    while(1)
#else
    Debug(libcw_do.off());
    for (math::CubicPolynomial& cubic : cubics)
#endif
    {
#if INTERACTIVE
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      double C0_s = slider_c0.value();
      math::CubicPolynomial cubic(C0_s, -3, 0, 1);

      plot::BezierFitter plot_bezier_fitter;
#endif
      int iterations;
      double initial_guess;
#if INTERACTIVE
      int n = get_roots(plot, second_layer, curve_line_style, &plot_bezier_fitter, cubic, roots, initial_guess, iterations);
#else
      int n = cubic.get_roots(roots, iterations);
#endif
      total_iterations += iterations;
      ++count;
      std::cout << "iterations = " << iterations << "; " << (total_iterations / count) << std::endl;

#if INTERACTIVE
      std::array<double, 3> real_roots;
      math::CubicPolynomial cubic_real(C0_s, -3, 0, 1);
      cubic_real.get_roots(real_roots, iterations);
      Dout(dc::notice, "Real root found: " << real_roots[0] << "; initial guess: " << initial_guess);

      // Draw a vertical line at the first found root.
      auto plot_root_line = plot.create_line(second_layer, line_style({.line_color = color::green}), Point{roots[0], 0.0}, Direction::up);

      // Draw a line at the guess.
      auto plot_guess_line = plot.create_line(second_layer, line_style, Point{initial_guess, 0.0}, Direction::up);
      auto plot_guess2_line = plot2.create_line(second_layer2, line_style, Point{slider_c0.value(), 0.0}, Direction::up);
      auto plot_guess3_line = plot2.create_line(second_layer2, line_style, Point{0.0, diff_lambda(slider_c0.value()).y()}, Direction::right);
#endif

#if 0
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
#endif

#if INTERACTIVE
      Debug(libcw_do.off());
      plot::BezierFitter plot_curve;
      plot_curve.solve(
          rel_err_lambda,
          plot2.viewport());
      plot2.add_bezier_fitter(second_layer2, curve_line_style, plot_curve);
      Debug(libcw_do.on());

      Debug(libcw_do.off());
      plot::BezierFitter plot_curve3;
      plot_curve3.solve(
          SC0_lambda,
          plot3.viewport());
      plot3.add_bezier_fitter(second_layer3, curve_line_style, plot_curve3);
      Debug(libcw_do.on());

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // Block until a redraw is necessary (for example because the user moved a draggable object,
      // or wants to print the current drawing) then go to the top of loop for a redraw.
      if (!window.handle_input_events())
        break;          // Program must be terminated.
#endif
    }
#if !INTERACTIVE
    Debug(libcw_do.on());
    ASSERT(count == number_of_cubics);
    Dout(dc::notice, "Average number of iterations per cubic: " << (total_iterations / number_of_cubics));
#endif

#if INTERACTIVE
    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }
#endif

  Dout(dc::notice, "Leaving main()");
}
