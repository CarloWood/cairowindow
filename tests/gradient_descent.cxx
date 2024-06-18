#include "sys.h"
#include "Polynomial.h"
#include "QuadraticPolynomial.h"
#include "Sample.h"
#include "Scale.h"
#include "Approximation.h"
#include "History.h"
#include "LocalExtreme.h"
#include "KineticEnergy.h"
#include "HorizontalDirection.h"
#include "VerticalDirection.h"
#include "Algorithm.h"
#include "PlotSample.h"
#include "PlotHistory.h"
#include "PlotLocalExtreme.h"
#include "PlotKineticEnergy.h"
#include "PlotParabolaScale.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/symbolic/symbolic.h"
#include "utils/square.h"
#include "utils/AIAlert.h"
#include <sstream>
#include <thread>
#include <vector>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/debug_ostream_operators.h"
#endif

#define USE_SLIDERS 0

using utils::has_print_on::operator<<;

class FunctionBase
{
 protected:
  symbolic::Symbol const& w_ = symbolic::Symbol::realize("w");
  symbolic::Symbol const& a_ = symbolic::Symbol::realize("a");
  symbolic::Symbol const& b_ = symbolic::Symbol::realize("b");
  symbolic::Symbol const& c_ = symbolic::Symbol::realize("c");
  symbolic::Symbol const& d_ = symbolic::Symbol::realize("d");
  symbolic::Symbol const& e_ = symbolic::Symbol::realize("e");
};

#if 0
class Function : FunctionBase
{
 public:
  static constexpr double w_0 = 1.0;
  static constexpr double w_min = -3.0;
  static constexpr double w_max = 6.0;

 private:
  using Expression = symbolic::Expression;

  Expression const& function_ = a_ + b_ * w_ + c_ * (w_^2) + d_ * (w_^3) + e_ * (w_^4);
  Expression const& derivative_ = function_.derivative(w_);

 public:
  Function()
  {
    a_ = 25.0;
    b_ = 25.0;
    c_ = -10.0;
    d_ = -10.0;
    e_ = 2.0;
  }

 public:
  std::string as_string() const
  {
    std::ostringstream oss;
    function_.print_on(oss);
    return oss.str();
  }

  double operator()(double w) const
  {
    w_ = w;
    return function_.evaluate();
  }

  double derivative(double w) const
  {
    w_ = w;
    return derivative_.evaluate();
  }
};
#elif 1
class Function : FunctionBase
{
 public:
  static constexpr double w_0 = -51.0;
  static constexpr double w_min = -60.0; //-80.0;
  static constexpr double w_max = -40.0; //20.0;

 private:
  using Expression = symbolic::Expression;
  using Constant = symbolic::Constant;
  using Symbol = symbolic::Symbol;
  using Func = symbolic::Function;

  static Constant const tp;
  Symbol const& amplitude_ = symbolic::Symbol::realize("amplitude");
  Symbol const& level_ = symbolic::Symbol::realize("level");
  Symbol const& phase_ = symbolic::Symbol::realize("phase");

 private:
  Func const& sigmoid_ = Func::realize("sigmoid", exp(3 * (w_ + tp)) / (1 + exp(3 * (w_ + tp))));

  Expression const& function_ = (1.0 - sigmoid_) * (a_ + b_ * w_ + c_ * (w_^2)) +
    (sigmoid_ * (amplitude_ * exp((tp - w_) / 10) * sin(d_ * w_ + phase_) + level_));

  Expression const& derivative_ = function_.derivative(w_);

 public:
  Function()
  {
    a_ = 14.6;
    b_ = 3.15;
    c_ = 0.451;
    d_ = 0.5;
    amplitude_ = 0.012027;
    level_ = 1878.38;
    phase_ = 1.91892;
  }

  std::string as_string() const
  {
    std::ostringstream oss;
    function_.print_on(oss);
    return oss.str();
  }

  double operator()(double w) const
  {
    w_ = w;
    sigmoid_.reset_evaluation();
    return function_.evaluate();
  }

  double derivative(double w) const
  {
    w_ = w;
    sigmoid_.reset_evaluation();
    return derivative_.evaluate();
  }

#if USE_SLIDERS
 public:
  void set_sliders(double amplitude, double level, double phase)
  {
    amplitude_ = amplitude;
    level_ = level;
    phase_ = phase;
    sigmoid_.reset_evaluation();
  }
#endif
};

//static
symbolic::Constant const Function::tp = symbolic::Constant::realize(55);

#elif 0
class Function : FunctionBase
{
 public:
  static constexpr double w_0 = 12.0;
  static constexpr double w_min = -20.0;
  static constexpr double w_max = 80.0;

 private:
  //symbolic::Expression const& function_ = a_ + b_ * w_ + c_ * (w_^2) + d_ * (w_^3) + e_ * (w_^4);
  symbolic::Expression const& function_ = symbolic::sin(a_ + b_ * w_) + d_ * ((w_ - c_)^2);
  symbolic::Expression const& derivative_ = function_.derivative(w_);

 public:
  Function()
  {
    //    4        3
    //  w    130⋅w          2
    //  ── - ────── + 2200⋅w  - 32000⋅w + 10000
    //  4      3

    a_ = 15.0;
    b_ = 0.3;
    c_ = 49.4;
    d_ = 0.0001;
#if 0
    a_ = 10000.0;
    b_ = -32000.0;
    c_ = 2200.0;
    d_ = -130.0 / 3.0;
    e_ = 0.25;
#endif
  }

  std::string as_string() const
  {
    std::ostringstream oss;
    function_.print_on(oss);
    return oss.str();
  }

  double operator()(double w) const
  {
    w_ = w;
    return function_.evaluate();
  }

  double derivative(double w) const
  {
    w_ = w;
    return derivative_.evaluate();
  }
};
#else
class Function
{
 public:
  static constexpr double w_0 = 12.0;
  static constexpr double w_min = -31.0;
  static constexpr double w_max = 16.0;

 private:
  using Expression = symbolic::Expression;
  using Constant = symbolic::Constant;
  using Symbol = symbolic::Symbol;
  using Func = symbolic::Function;

 private:
  Symbol const& w_ = Symbol::realize("w");
  Symbol const& a1_ = Symbol::realize("a1");
  Symbol const& a2_ = Symbol::realize("a2");
  Symbol const& b1_ = Symbol::realize("b1");
  Symbol const& b2_ = Symbol::realize("b2");
  Symbol const& c1_ = Symbol::realize("c1");
  Symbol const& c2_ = Symbol::realize("c2");

  Func const& sigmoid_ = Func::realize("sigmoid", exp(5 * w_) / (1 + exp(5 * w_)));

  Expression const& function_ = (a1_ + (a2_ - a1_) * sigmoid_ + (b1_ + (b2_ - b1_) * sigmoid_) * w_ + (c1_ + (c2_ - c1_) * sigmoid_) * (w_^2));
  Expression const& derivative_ = function_.derivative(w_);

 public:
  Function()
  {
    a1_ = 15.0;
    b1_ = 3.1;
    c1_ = 0.2;
    a2_ = 14.6;
    b2_ = 3.15;
    c2_ = 0.451;
  }

  std::string as_string() const
  {
    std::ostringstream oss;
    function_.print_on(oss);
    return oss.str();
  }

  double operator()(double w) const
  {
    w_ = w;
    sigmoid_.reset_evaluation();
    double result = function_.evaluate();
    return result;
  }

  double derivative(double w) const
  {
    w_ = w;
    sigmoid_.reset_evaluation();
    return derivative_.evaluate();
  }
};
#endif

//static
cairowindow::draw::ConnectorStyle const PlotParabolaScale::s_indicator_style{{.line_width = 1}};

//static
cairowindow::draw::ConnectorStyle const PlotKineticEnergy::s_indicator_style{{.line_width = 1}};

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  Function L;

  double const w_0 = Function::w_0;
  double const w_min = Function::w_min;
  double const w_max = Function::w_max;

  int const steps = 100;

  try
  {
    using namespace cairowindow;
    using namespace gradient_descent;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Gradient descent of " + L.as_string(), 2*600, 2*450);

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
//      L_min = -100;
//      L_max = 60.0;
    }

    // Create and draw plot area.
    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "(1−σ(3(w+55))) (14.6+3.15w+0.451w²) + σ(3(w+55))(0.012 exp((55-w)/10) sin(0.5w+1.92)+1878)", {},
        "w", {},
        "L", {});
    plot.set_xrange({w_min, w_max});
    plot.set_yrange({L_min, L_max});
    plot.add_to(background_layer, false);

    draw::PointStyle point_style({.color_index = 31, .filled_shape = 1});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::LineStyle derivative_line_style({.line_color = color::turquoise, .line_width = 1.0});

#if !USE_SLIDERS
    BezierFitter L_fitter([&L](double w) -> Point { return {w, L(w)}; }, plot.viewport());
    auto plot_curve = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(L_fitter));
#endif

    gradient_descent::Algorithm gda(0.1, L_max, plot, second_layer, point_style, label_style);

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

    // Loop over iterations of w.
    for (Weight w(w_0);;)
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
        break;
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

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
