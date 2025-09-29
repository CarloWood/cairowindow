#include "sys.h"
#include "AlgorithmEvent.h"
#include "math/Polynomial.h"
#include "math/QuadraticPolynomial.h"
#include "gradient_descent/Sample.h"
#include "gradient_descent/Scale.h"
#include "gradient_descent/Approximation.h"
#include "gradient_descent/History.h"
#include "gradient_descent/LocalExtreme.h"
#include "gradient_descent/KineticEnergy.h"
#include "gradient_descent/HorizontalDirection.h"
#include "gradient_descent/ExtremeType.h"
#include "gradient_descent/Algorithm.h"
#include "gradient_descent/History.h"
#include "gradient_descent/LocalExtreme.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/symbolic/symbolic.h"
#include "utils/square.h"
#include "utils/AIAlert.h"
#include <sstream>
#include <thread>
#include <vector>
#include <random>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/debug_ostream_operators.h"
#endif

#define USE_SLIDERS 0

using utils::has_print_on::operator<<;

class Function
{
 private:
  double w_0_;
  double w_min_;
  double w_max_;
  math::CubicPolynomial<double> cubic_;

 public:
  Function(std::mt19937& engine)
  {
    std::uniform_real_distribution<double> edge_dist(-20.0, 20.0);
    w_min_ = edge_dist(engine);
    w_max_ = edge_dist(engine);
    if (w_max_ < w_min_)
      std::swap(w_min_, w_max_);

    std::uniform_real_distribution<double> L_dist(-20.0, 20.0);
    double Lle = L_dist(engine);
    double Lre = L_dist(engine);

    std::uniform_real_distribution<double> slope_dist(-0.4 * M_PI, 0.4 * M_PI);
    double la = slope_dist(engine);
    double ra = slope_dist(engine);
    double dle = std::atan(la);
    double dre = std::atan(ra);

    cubic_.initialize(w_min_, Lle, dle, w_max_, Lre, dre);

    std::array<double, 2> extrema;
    int number_of_extrema = cubic_.get_extrema(extrema);

    if (number_of_extrema == 2)
    {
      double width = std::max(w_max_, extrema[1]) - std::min(w_min_, extrema[0]);
      w_min_ = std::min(w_min_, extrema[0] - 0.1 * width);
      w_max_ = std::max(w_max_, extrema[1] + 0.1 * width);
    }

    std::uniform_real_distribution<double> w0_dist(w_min_, w_max_);
    w_0_ = w0_dist(engine);
  }

  double operator()(double w) const
  {
    return cubic_(w);
  }

  double derivative(double w) const
  {
    return cubic_.derivative(w);
  }

  std::string as_string() const
  {
    return "Generated function";
  }

  double w_0() const { return w_0_; }
  double w_min() const { return w_min_; }
  double w_max() const { return w_max_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << as_string();
  }
#endif
};

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  int seed = 10248;
  Dout(dc::notice, "seed = " << seed);
  std::mt19937 engine(seed);

  try
  {
    using namespace cairowindow;
    using namespace gradient_descent;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("gradient_descent", 2*600, 2*450);

    // Create a new layer with a white background.
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

    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::LineStyle derivative_line_style({.line_color = color::turquoise, .line_width = 1.0});

#if USE_SLIDERS
    // amplitude = 0.012027, level = 1878.38, phase = 1.91892
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    auto slider_amplitude = plot.create_slider(second_layer, {978, 83, 7, 400}, 0.012027, 0.01, 0.04);
    auto slider_amplitude_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "amplitude");
    auto slider_level = plot.create_slider(second_layer, {1028, 83, 7, 400}, 1878.38, 1500, 2500);
    auto slider_level_label = plot.create_text(second_layer, slider_style, Pixel{1028, 483}, "level");
    auto slider_phase = plot.create_slider(second_layer, {1078, 83, 7, 400}, 1.91892, 0.0, 2.0 * M_PI);
    auto slider_phase_label = plot.create_text(second_layer, slider_style, Pixel{1078, 483}, "phase");
#endif

    for (int t = 0; t < 100; ++t)
    {
      Function L(engine);

      double const w_0 = L.w_0();
      double const w_min = L.w_min();
      double const w_max = L.w_max();

      int const steps = 100;
      double L_min = L(w_min);
      double L_max = L_min;
      {
        double w = w_min;
        double delta_w = (w_max - w_min) / (steps - 1);
        for (int i = 0; i < steps; ++i)
        {
          double val = L(w);
          L_min= std::min(L_min, val);
          L_max= std::max(L_max, val);
          w += delta_w;
        }
      }

      // Create and draw plot area.
      plot::Plot plot(window.geometry(), draw::PlotAreaStyle({.color = color::orange}),
          L.as_string(), {},
          "w", {},
          "L", {});
      plot.set_xrange({w_min, w_max});
      plot.set_yrange({L_min, L_max});
      plot.add_to(background_layer, false);

#if !USE_SLIDERS
      BezierFitter L_fitter([&L](double w) -> Point { return {w, L(w)}; }, plot.viewport());
      auto plot_curve = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(L_fitter));
#endif

      gradient_descent::Algorithm gda(0.1, L_max);

#ifdef CWDEBUG
      AlgorithmEvent algorithm_event(plot, second_layer);
      auto algorithm_event_handle = gda.event_server().request(algorithm_event, &AlgorithmEvent::callback);
#endif

      // Loop over iterations of w.
      bool next_curve = false;
      for (Weight w(w_0); !next_curve;)
      {
        // Suppress immediate updating of the window for each created item, in order to avoid flickering.
        window.set_send_expose_events(false);

#if USE_SLIDERS
        L.set_sliders(slider_amplitude.value(), slider_level.value(), slider_phase.value());
        Dout(dc::notice,
            "amplitude = " << slider_amplitude.value() << ", level = " << slider_level.value() << ", phase = " << slider_phase.value());

        BezierFitter L_fitter([&L](double w) -> Point { return {w, L(w)}; }, plot.viewport());
        auto plot_curve = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(L_fitter));
#else
        // Replace w with a new point until the global minimum has been reached.
        if (!gda(w, L(w), L.derivative(w)))
          next_curve = true;
#endif

        // Flush all expose events related to the drawing done above.
        window.set_send_expose_events(true);

        // Block until a key is pressed.
        if (!window.handle_input_events())
          break;          // Program must be terminated.

        Dout(dc::notice, "------------------------------------");
      }

#if !USE_SLIDERS
      if (gda.success())
        Dout(dc::notice, "Found global minimum " << gda.minimum());
#endif

      Debug(algorithm_event_handle.cancel());
    }

    window.close();
    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
