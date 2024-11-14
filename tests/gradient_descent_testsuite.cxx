#include "sys.h"
#include "AlgorithmEvent.h"
#include "gradient_descent/Algorithm.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Range.h"
#include "cairowindow/Plot.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/symbolic/symbolic.h"
#include "utils/almost_equal.h"
#include <thread>
#include "debug.h"

using namespace cairowindow;

struct EnableDrawing;

struct Function
{
  symbolic::Function const& function_;
  symbolic::Symbol const& symbol_;
  std::vector<symbolic::Function const*> deps_;

  template<typename... Deps>
  Function(symbolic::Symbol const& symbol, symbolic::Function const& function, Deps&&... deps) :
    function_(function), symbol_(symbol), deps_({&deps...}) { }

  double operator()(double x) const
  {
    symbol_ = x;
    reset_evaluation();
    return function_.evaluate();
  }

  double derivative(double x) const
  {
    symbol_ = x;
    reset_evaluation();
    return function_.derivative(symbol_).evaluate();
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    function_.print_on(os);
  }
#endif

 private:
  friend struct EnableDrawing;
  void reset_evaluation() const
  {
    function_.reset_evaluation();
    for (symbolic::Function const* dep : deps_)
      dep->reset_evaluation();
  }

#ifdef CWDEBUG
  std::string to_string() const
  {
    return function_.to_string();
  }
#endif
};

class EnableDrawing
{
 public:
  static constexpr draw::LineStyle curve_line_style{{.line_width = 1.0}};

 private:
  cairowindow::Window window;
  boost::intrusive_ptr<Layer> background_layer;
  boost::intrusive_ptr<Layer> second_layer;
  std::thread event_loop;
  Range L_min_max;
  plot::Plot plot;
  plot::BezierFitter plot_curve;
#ifdef CWDEBUG
  AlgorithmEvent algorithm_event;
  events::RequestHandle<gradient_descent::AlgorithmEventType> algorithm_event_handle;
#endif

 public:
  EnableDrawing(gradient_descent::Algorithm* algorithm, Function const& L, double w_min, double w_max);
  ~EnableDrawing();

  static Range get_L_min_max(Function const& L, double w_min, double w_max)
  {
    double L_min = L(w_min);
    double L_max = L_min;
    {
      int const steps = 100;
      double w = w_min;
      double delta_w = (w_max - w_min) / (steps - 1);
      for (int i = 0; i < steps; ++i)
      {
        double Lw = L(w);
        L_min= std::min(L_min, Lw);
        L_max= std::max(L_max, Lw);
        w += delta_w;
      }
    }
    return {L_min, L_max};
  }

  void wait()
  {
    // Block until a key is pressed.
    if (!window.handle_input_events())
      Dout(dc::notice, "Algorithm Terminated");
  }
};

EnableDrawing::EnableDrawing(gradient_descent::Algorithm* algorithm, Function const& L, double w_min, double w_max) :
  window("Gradient descent of "
#ifdef CWDEBUG
      + L.to_string(),
#else
      "L",
#endif
      2*600, 2*450),
  background_layer(window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"))),
  second_layer(window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"))),
  event_loop([&](){
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    }),
  L_min_max(get_L_min_max(L, w_min, w_max)),
  plot(window.geometry(), { .grid = {.color = color::orange} },
#ifdef CWDEBUG
        L.to_string(),
#else
        "L",
#endif
        {}, "w", {}, "L", {})
        COMMA_DEBUG_ONLY(algorithm_event(plot, second_layer),
            algorithm_event_handle(algorithm->event_server().request(algorithm_event, &AlgorithmEvent::callback)))
{
  plot.set_xrange({w_min, w_max});
  plot.set_yrange(L_min_max);
  plot.add_to(background_layer, false);

  plot_curve.solve([&L](double w) -> Point { return {w, L(w)}; }, plot.viewport());
  plot.add_bezier_fitter(second_layer, curve_line_style, plot_curve);

  L.reset_evaluation();
}

EnableDrawing::~EnableDrawing()
{
#ifdef CWDEBUG
  algorithm_event_handle.cancel();
#endif
  window.close();
  event_loop.join();
}

class Algorithm : public gradient_descent::Algorithm
{
 private:
  std::unique_ptr<EnableDrawing> enable_drawing_;

 public:
  using gradient_descent::Algorithm::Algorithm;

  void enable_drawing(Function const& L, double w_min, double w_max)
  {
    enable_drawing_ = std::make_unique<EnableDrawing>(this, L, w_min, w_max);
  }

  bool operator()(gradient_descent::Weight& w, double Lw, double dLdw)
  {
    bool result = gradient_descent::Algorithm::operator()(w, Lw, dLdw);
    if (enable_drawing_)
      enable_drawing_->wait();
    return result;
  }
};

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  using Weight = gradient_descent::Weight;
  using Scale = gradient_descent::Scale;
  using HorizontalDirection = gradient_descent::HorizontalDirection;
  using Sample = gradient_descent::Sample;
  using Approximation = gradient_descent::Approximation;
  using ExtremeType = gradient_descent::ExtremeType;

  // Default values.
  constexpr double L_max = 100.0;
  constexpr double learning_rate = 0.125;

#if 0
  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of zero ***");
  {
    // Algorithm: [one sample, derivative is zero, hdirection is unknown]
    constexpr double w0 = 13.0;

    Algorithm gda(learning_rate, L_max);
    // Drop in with a zero derivative.
    Weight w = w0;
    gda(w, 50.0, 0.0);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero, hdirection is unknown");

    // EXPECTED: learning_rate was added.
    ASSERT(w == w0 + learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of almost zero (less than one million-th w0) ***");
  {
    // Algorithm: [one sample, derivative is zero, hdirection is unknown]
    constexpr double w0 = 13.0;

    Algorithm gda(learning_rate, L_max);
    // Drop in with an almost zero derivative.
    Weight w = w0;
    gda(w, 50.0, w0 * 0.999e-6);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero, hdirection is unknown");

    // EXPECTED: learning_rate was added.
    ASSERT(w == w0 + learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of almost zero (less than epsilon), but w0 is also zero ***");
  {
    // Algorithm: [one sample, derivative is zero, hdirection is unknown]
    constexpr double w0 = 0.0;

    Algorithm gda(learning_rate, L_max);
    // Drop in with an almost zero derivative.
    Weight w = w0;
    gda(w, 50.0, 0.999 * Scale::epsilon);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero, hdirection is unknown");

    // EXPECTED: learning_rate was added.
    ASSERT(w == w0 + learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a small but non-zero derivative while w0 is zero ***");
  {
    // Algorithm: [one sample, gradient descent]
    constexpr double w0 = 0.0;
    constexpr double dLdw = 1.5e-8;

    Algorithm gda(learning_rate, L_max);
    // Drop in with an almost zero derivative.
    Weight w = w0;
    gda(w, 50.0, dLdw);
    ASSERT(gda.algorithm_str() == "one sample, gradient descent");

    // EXPECTED: learning_rate * derivative was subtracted.
    ASSERT(w == w0 - learning_rate * dLdw);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of zero (while hdirection is known) ***");
  std::array<HorizontalDirection, 2> hdirections = {
    HorizontalDirection::left,
    HorizontalDirection::right
  };
  for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [one sample, derivative is zero]
    constexpr double w0 = 13.0;
    HorizontalDirection const hdirection = hdirections[hdi];

    Dout(dc::notice, "* hdirection = " << hdirection);

    Algorithm gda(learning_rate, L_max);
    ASSERT(gda.debug_small_step() == 0.0);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, gda.debug_next_extreme_type(), 0.0);

    // Drop in with a zero derivative.
    Weight w = w0;
    gda(w, 50.0, 0.0);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero");

    // EXPECTED: hdirection * learning_rate was added.
    ASSERT(gda.debug_hdirection() == hdirection);
    ASSERT(w == w0 + static_cast<int>(hdirection) * learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of almost zero (less than one million-th w0) (while hdirection is known) ***");
  for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [one sample, derivative is zero]
    constexpr double w0 = 13.0;
    HorizontalDirection const hdirection = hdirections[hdi];

    Dout(dc::notice, "* hdirection = " << hdirection);

    Algorithm gda(learning_rate, L_max);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, gda.debug_next_extreme_type(), 0.0);

    // Drop in with an almost zero derivative.
    Weight w = w0;
    gda(w, 50.0, w0 * 0.999e-6);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero");

    // EXPECTED: hdirection * learning_rate was added.
    ASSERT(gda.debug_hdirection() == hdirection);
    ASSERT(w == w0 + static_cast<int>(hdirection) * learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a derivative of almost zero (less than epsilon), but w0 is also zero (while hdirection is known) ***");
  for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [one sample, derivative is zero]
    constexpr double w0 = 0.0;
    HorizontalDirection const hdirection = hdirections[hdi];

    Dout(dc::notice, "* hdirection = " << hdirection);

    Algorithm gda(learning_rate, L_max);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, gda.debug_next_extreme_type(), 0.0);

    // Drop in with an almost zero derivative.
    Weight w = w0;
    gda(w, 50.0, 0.999 * Scale::epsilon);
    ASSERT(gda.algorithm_str() == "one sample, derivative is zero");

    // EXPECTED: hdirection * learning_rate was added.
    ASSERT(gda.debug_hdirection() == hdirection);
    ASSERT(w == w0 + static_cast<int>(hdirection) * learning_rate);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with a small but non-zero derivative while w0 is zero (while hdirection is known) ***");
  for (int sign_of_dLdw = -1; sign_of_dLdw <= 1; sign_of_dLdw +=2)
    for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [one sample, same direction]
    constexpr double w0 = 0.0;
    HorizontalDirection const hdirection = hdirections[hdi];
    double const dLdw = sign_of_dLdw * 1.5e-8;

    Dout(dc::notice, "* hdirection = " << hdirection << "; dLdw = " << dLdw);

    Algorithm gda(learning_rate, L_max);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, gda.debug_next_extreme_type(), 0.0);

    // Drop in with an almost zero derivative.
    Weight w = w0;
    gda(w, 50.0, dLdw);
    ASSERT(gda.algorithm_str() == "one sample, same direction");

    // EXPECTED: hdirection_ * abs(learning_rate_ * dLdW) was added.
    ASSERT(gda.debug_hdirection() == hdirection);
    ASSERT(w == w0 + static_cast<int>(hdirection) * std::abs(learning_rate * dLdw));
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: starting with zero derivative but known small_step (and thus while hdirection is known) ***");
  for (int hdi = 0; hdi < hdirections.size(); ++hdi)
  {
    // Algorithm: [small step]
    constexpr double w0 = 13.0;
    constexpr double small_step = 0.12345;
    HorizontalDirection const hdirection = hdirections[hdi];

    Dout(dc::notice, "* hdirection = " << hdirection);

    Algorithm gda(learning_rate, L_max);
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, ExtremeType::maximum, small_step);

    // Drop in with a zero derivative.
    Weight w = w0;
    gda(w, 50.0, 0.0);
    ASSERT(gda.algorithm_str() == "small step");

    // EXPECTED: small_step was subtracted (rather randomly, because we are looking for a maximum).
    ASSERT(utils::almost_equal(static_cast<double>(w), w0 - small_step, 10e-15));
  }
#endif

#if 0
  //==========================================================================
  Dout(dc::notice, "*** TEST: jump to minimum ***");
  {
    constexpr double w0 = 13.0;

    Algorithm gda(learning_rate, L_max);
    Weight w = w0;
    symbolic::Symbol const& x = symbolic::Symbol::realize("w");
    symbolic::Constant const& c4 = symbolic::Constant::realize(1, 10000000);
    symbolic::Function const& sL = symbolic::Function::realize("L", (2 - 5 * x + (x^2)) - c4 * (x^4));
    constexpr double cp_w = 2.5001550590858685;
    Function L(x, sL);

    //gda.enable_drawing(L, -20.0, 20.0);

    Dout(dc::notice, "L = " << fulldef << L);

    double dLdw = L.derivative(w);
    gda(w, L(w), dLdw);                                 // 10.3750010985 [one sample, gradient descent]
    ASSERT(gda.algorithm_str() == "one sample, gradient descent");
    ASSERT(utils::almost_equal(static_cast<double>(w), w0 - learning_rate * dLdw, 1e-15));
    gda(w, L(w), L.derivative(w));                      // 2.50015505909 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "initial find_extreme jump");
    // This is the vertex.
    double const vertex = w;
    ASSERT(utils::almost_equal(vertex, cp_w, 1e-15));
    dLdw = L.derivative(w);
    ASSERT(dLdw < 1e-3);
    gda(w, L(w), 0.0);                                 // -7.99968988182 [keep going (no zeroes)]
    ASSERT(gda.algorithm_str() == "keep going (no zeroes)");
    double past_extreme = w;
    // The same value should not be returned twice (this requires a special-case check).
    ASSERT(!utils::almost_equal(past_extreme, cp_w, 1e-3));
    // Instead, we should have went further to the left.
    ASSERT(past_extreme < -6.0);
    bool aborted;
    do
    {
      aborted = !gda(w, L(w), L.derivative(w));         // -6241.8483809533418 [small step]
                                                        // -13752.619901104652 [continue same direction]
                                                        // 13 [best minimum, opposite direction]
    }
    while (w < past_extreme);
    if (!aborted)
    {
      ASSERT(gda.algorithm_str() == "best minimum, opposite direction");
      // Since this is the sample that it the furthest away from the critical point and clearly
      // part of the cubic approximation, that is where we expect to end up again
      // when adding the scale to the (best) minimum.
      ASSERT(utils::almost_equal(static_cast<double>(w), static_cast<double>(w0), 10e-9));

      // At this point, small_step must be zero and hdirection right.
      ASSERT(gda.debug_hdirection() == HorizontalDirection::right);
      ASSERT(utils::almost_equal(gda.debug_small_step(), w0 - vertex, 10e-6));
    }
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: connected parabolas ***");
  {
    constexpr double w0 = 12.0;
    constexpr double learning_rate = 0.1;
    constexpr double L_max = 180.456;

    Algorithm gda(learning_rate, L_max);
    Weight w = w0;

    symbolic::Constant const& two = symbolic::Constant::realize(2);
    symbolic::Symbol const& x = symbolic::Symbol::realize("w");
    symbolic::Function const& sigmoid = symbolic::Function::realize("sigmoid", exp(5 * x) / (1 + exp(5 * x)));
    symbolic::Constant const& a1 = symbolic::Constant::realize(15);
    symbolic::Constant const& b1 = symbolic::Constant::realize(31, 10);
    symbolic::Constant const& c1 = symbolic::Constant::realize(2, 10);
    double vertex1 = (-b1 / (two * c1)).evaluate();
    symbolic::Constant const& a2 = symbolic::Constant::realize(146, 10);
    symbolic::Constant const& b2 = symbolic::Constant::realize(315, 100);
    symbolic::Constant const& c2 = symbolic::Constant::realize(451, 1000);
    double vertex2 = (-b2 / (two * c2)).evaluate();
    symbolic::Function const& sL = symbolic::Function::realize("L",
        a1 + (a2 - a1) * sigmoid + (b1 + (b2 - b1) * sigmoid) * x + (c1 + (c2 - c1) * sigmoid) * (x^2));
    Function L(x, sL, sigmoid);

//    gda.enable_drawing(L, -35.0, 20.0);

    double dLdw = L.derivative(w);
    gda(w, L(w), dLdw);
    ASSERT(gda.algorithm_str() == "one sample, gradient descent");
    ASSERT(utils::almost_equal(static_cast<double>(w), w0 - learning_rate * dLdw, 1e-15));
    gda(w, L(w), L.derivative(w));                      // -3.49223946783 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "initial find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), vertex2, 1e-9));
    gda(w, L(w), L.derivative(w));                      // -7.65806821678 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), -7.65806821678, 1e-9));
    gda(w, L(w), L.derivative(w));                      // -7.75000002673 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), vertex1, 1e-8));
    gda(w, L(w), L.derivative(w));                      // -12.0077604322 [keep going (no zeroes)]
    ASSERT(gda.algorithm_str() == "keep going (no zeroes)");
    ASSERT(utils::almost_equal(static_cast<double>(w), -12.0077604322, 1e-8));
    gda(w, L(w), L.derivative(w));                      // -16.2655209145 [keep going]
    ASSERT(gda.algorithm_str() == "keep going");
    ASSERT(utils::almost_equal(static_cast<double>(w), -16.2655209145, 1e-8));
    gda(w, L(w), L.derivative(w));                      // -23.3617883757 [keep going]
    ASSERT(gda.algorithm_str() == "keep going");
    ASSERT(utils::almost_equal(static_cast<double>(w), -23.3617883757, 1e-8));
    gda(w, L(w), L.derivative(w));                      // -36.5405707968 [keep going]
    ASSERT(gda.algorithm_str() == "keep going");
    ASSERT(utils::almost_equal(static_cast<double>(w), -36.5405707968, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 5.42878239433 [best minimum, opposite direction]
    ASSERT(gda.algorithm_str() == "best minimum, opposite direction");
    ASSERT(utils::almost_equal(static_cast<double>(w), 5.42878239433, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 18.6075648154 [keep going (bad extreme)]
    ASSERT(gda.algorithm_str() == "keep going (bad extreme)");
    ASSERT(utils::almost_equal(static_cast<double>(w), 18.6075648154, 1e-8));
    bool finished = !gda(w, L(w), L.derivative(w));
    ASSERT(finished);
    ASSERT(gda.success());
    Sample const& result = gda.minimum();
    ASSERT(utils::almost_equal(result.w(), vertex1, 1e-8));
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: parabola plus sin ***");
  {
    constexpr double w0 = 12.0;
    constexpr double learning_rate = 0.1;
    constexpr double L_max = 0;

    Algorithm gda(learning_rate, L_max);
    Weight w = w0;

    symbolic::Symbol const& x = symbolic::Symbol::realize("w");
    symbolic::Constant const& a = symbolic::Constant::realize(15);
    symbolic::Constant const& b = symbolic::Constant::realize(3, 10);
    symbolic::Constant const& c = symbolic::Constant::realize(494, 10);
    symbolic::Constant const& d = symbolic::Constant::realize(1, 10000);
    symbolic::Function const& sL = symbolic::Function::realize("L", sin(a + b * x) + d * ((x - c)^2));
    Function L(x, sL);

//    gda.enable_drawing(L, -20.0, 80.0);

    gda(w, L(w), L.derivative(w));                      // 11.9716773342 [one sample, gradient descent]
    ASSERT(gda.algorithm_str() == "one sample, gradient descent");
    ASSERT(utils::almost_equal(static_cast<double>(w), 11.9716773342, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 8.12417155659 [initial find_extreme jump]
    ASSERT(gda.algorithm_str() == "initial find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), 8.12417155659, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 7.73480664643 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), 7.73480664643, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 7.68864155893 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), 7.68864155893, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 15.7327395835 [best zero]
    ASSERT(gda.algorithm_str() == "best zero");
    ASSERT(utils::almost_equal(static_cast<double>(w), 15.7327395835, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 17.667571638 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), 17.667571638, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 17.9850634722 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), 17.9850634722, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 17.9980425743 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), 17.9980425743, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 26.2673319378 [best zero]
    ASSERT(gda.algorithm_str() == "best zero");
    ASSERT(utils::almost_equal(static_cast<double>(w), 26.2673319378, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 29.0145983989 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), 29.0145983989, 1e-8));
    gda(w, L(w), L.derivative(w));                      // 28.59944683 [find_extreme jump]
    ASSERT(gda.algorithm_str() == "find_extreme jump");
    ASSERT(utils::almost_equal(static_cast<double>(w), 28.59944683, 1e-8));
#if 0
    bool finished = !gda(w, L(w), L.derivative(w));
    ASSERT(finished);
    ASSERT(gda.success());
    Sample const& result = gda.minimum();
    ASSERT(utils::almost_equal(result.w(), vertex1, 1e-8));
#endif
  }
#endif

  //==========================================================================
  Dout(dc::notice, "*** TEST: parabola connected to dampened sin ***");
  {
    constexpr double w0 = -10.0;
    constexpr double learning_rate = 0.1;
    constexpr double L_max = 2649;

    Algorithm gda(learning_rate, L_max);
    Weight w = w0;

    symbolic::Symbol const& x = symbolic::Symbol::realize("w");
    symbolic::Constant const& tp = symbolic::Constant::realize(55);
    symbolic::Function const& sigmoid = symbolic::Function::realize("sigmoid", exp(3 * (x + tp)) / (1 + exp(3 * (x + tp))));
    symbolic::Constant const& a = symbolic::Constant::realize(146, 10);
    symbolic::Constant const& b = symbolic::Constant::realize(315, 100);
    symbolic::Constant const& c = symbolic::Constant::realize(451, 1000);
    symbolic::Constant const& d = symbolic::Constant::realize(5, 10);
    symbolic::Constant const& amplitude = symbolic::Constant::realize(12027, 1000000);
    symbolic::Constant const& level = symbolic::Constant::realize(187838, 100);
    symbolic::Constant const& phase = symbolic::Constant::realize(191892, 100000);
    symbolic::Function const& sL = symbolic::Function::realize("L",
        (1.0 - sigmoid) * (a + b * x + c * (x^2)) + (sigmoid * (amplitude * exp((tp - x) / 10) * sin(d * x + phase) + level)));
    Function L(x, sL, sigmoid);

    //gda.enable_drawing(L, -52.0, -43.0);
    gda.enable_drawing(L, -80.0, 20.0);

// FIXME: Currently core dumps
    while (gda(w, L(w), L.derivative(w)))
    {
      Dout(dc::notice, "-------------------------------------------");
    }

    ASSERT(gda.success());
    Sample const& result = gda.minimum();
    //ASSERT(utils::almost_equal(result.w(), vertex1, 1e-8));
    Dout(dc::notice, "Global minimum: " << result);
  }

#if 0
  //==========================================================================
  Dout(dc::notice, "*** TEST: getting extra samples ***");
  {
    constexpr double w0 = 1.0;
    constexpr double learning_rate = 0.007557;
    constexpr double L_max = 1328.89;

    Algorithm gda(learning_rate, L_max);
    Weight w = w0;

    symbolic::Symbol const& x = symbolic::Symbol::realize("w");
    symbolic::Constant const& a = symbolic::Constant::realize(25);
    symbolic::Constant const& b = symbolic::Constant::realize(25);
    symbolic::Constant const& c = symbolic::Constant::realize(-10);
    symbolic::Constant const& d = symbolic::Constant::realize(-10);
    symbolic::Constant const& e = symbolic::Constant::realize(2);
    symbolic::Function const& sL = symbolic::Function::realize("L", a + b * x + c * (x^2) + d * (x^3) + e * (x^4));
    Function L(x, sL);

    gda.enable_drawing(L, -3.0, 7.2);

    while (gda(w, L(w), L.derivative(w)))
      ;
  }
#endif

  Dout(dc::notice, "Success!");
}
