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

  Function(symbolic::Function const& function, symbolic::Symbol const& symbol) :
    function_(function), symbol_(symbol) { }

  double operator()(double x) const
  {
    symbol_ = x;
    function_.reset_evaluation();
    return function_.evaluate();
  }

  double derivative(double x) const
  {
    symbol_ = x;
    function_.reset_evaluation();
    return function_.derivative(symbol_).evaluate();
  }

  void print_on(std::ostream& os) const
  {
    function_.print_on(os);
  }

 private:
  friend struct EnableDrawing;
  symbolic::Function const* operator->() const
  {
    return &function_;
  }
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
  AlgorithmEvent algorithm_event;
  events::RequestHandle<gradient_descent::AlgorithmEventType> algorithm_event_handle;
  plot::BezierFitter plot_curve;

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
  window("Gradient descent of " + L->to_string(), 2*600, 2*450),
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
        L->to_string(), {}, "w", {}, "L", {}),
  algorithm_event(plot, second_layer),
  algorithm_event_handle(algorithm->event_server().request(algorithm_event, &AlgorithmEvent::callback))
{
  plot.set_xrange({w_min, w_max});
  plot.set_yrange(L_min_max);
  plot.add_to(background_layer, false);

  plot_curve.solve([&L](double w) -> Point { return {w, L(w)}; }, plot.viewport());
  plot.add_bezier_fitter(second_layer, curve_line_style, plot_curve);

  L->reset_evaluation();
}

EnableDrawing::~EnableDrawing()
{
  algorithm_event_handle.cancel();
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

  // Default values.
  constexpr double L_max = 100.0;
  constexpr double learning_rate = 0.125;

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
    gda.debug_set_hdirection_next_extreme_type_small_step(hdirection, gda.debug_next_extreme_type(), small_step);

    // Drop in with a zero derivative.
    Weight w = w0;
    gda(w, 50.0, 0.0);
    ASSERT(gda.algorithm_str() == "small step");

    // EXPECTED: hdirection_ * small_step was added.
    ASSERT(gda.debug_hdirection() == hdirection);
    ASSERT(utils::almost_equal(static_cast<double>(w), w0 + static_cast<int>(hdirection) * small_step, 10e-15));
  }

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
    Function L(sL, x);

    gda.enable_drawing(L, -20.0, 20.0);

    Dout(dc::notice, "L = " << fulldef << L);

    double dLdw = L.derivative(w);
    gda(w, L(w), dLdw);                                 // 10.3750010985 [one sample, gradient descent]
    ASSERT(gda.algorithm_str() == "one sample, gradient descent");
    ASSERT(utils::almost_equal(static_cast<double>(w), w0 - learning_rate * dLdw, 1e-15));
    gda(w, L(w), L.derivative(w));                      // 2.5000015493572083 [jump to vertex]
    ASSERT(gda.algorithm_str() == "jump to vertex");
    // This is the vertex.
    double const vertex = w;
    ASSERT(utils::almost_equal(vertex, cp_w, 1e-15));
    dLdw = L.derivative(w);
    ASSERT(dLdw < 1e-3);
    gda(w, L(w), 0.0);                                 // -7.99968988182 [past extreme (no zeroes)]
    ASSERT(gda.algorithm_str() == "past extreme (no zeroes)");
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
      // Since this is the sample that it the furthest away from the vertex and clearly
      // part of the parabolic approximation, that is where we expect to end up again
      // when adding the scale to the (best) minimum.
      ASSERT(w == w0);

      // At this point, small_step must be zero and hdirection right.
      ASSERT(gda.debug_hdirection() == HorizontalDirection::right);
      ASSERT(utils::almost_equal(gda.debug_small_step(), w0 - vertex, 10e-6));
    }
  }

  Dout(dc::notice, "Success!");
}
