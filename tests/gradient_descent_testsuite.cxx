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

struct EnableDrawing
{
  static constexpr draw::LineStyle curve_line_style{{.line_width = 1.0}};

  cairowindow::Window window;
  boost::intrusive_ptr<Layer> background_layer;
  boost::intrusive_ptr<Layer> second_layer;
  std::thread event_loop;
  Range L_min_max;
  plot::Plot plot;
  AlgorithmEvent algorithm_event;
  events::RequestHandle<gradient_descent::AlgorithmEventType> algorithm_event_handle;
  plot::BezierFitter plot_curve;

  EnableDrawing(gradient_descent::Algorithm& algorithm, symbolic::Function const& L, symbolic::Symbol const& x,
      double w_min, double w_max);

  ~EnableDrawing();

  static Range get_L_min_max(symbolic::Function const& L, symbolic::Symbol const& x, double w_min, double w_max)
  {
    x = w_min;
    L.reset_evaluation();
    double L_min = L.evaluate();
    double L_max = L_min;
    {
      int const steps = 100;
      double w = w_min;
      double delta_w = (w_max - w_min) / (steps - 1);
      for (int i = 0; i < steps; ++i)
      {
        x = w;
        L.reset_evaluation();
        double val = L.evaluate();
        L_min= std::min(L_min, val);
        L_max= std::max(L_max, val);
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

EnableDrawing::EnableDrawing(gradient_descent::Algorithm& algorithm, symbolic::Function const& L, symbolic::Symbol const& x,
    double w_min, double w_max) :
  window("Gradient descent of " + L.to_string(), 2*600, 2*450),
  background_layer(window.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("background_layer"))),
  second_layer(window.create_layer<Layer>({} COMMA_DEBUG_ONLY("second_layer"))),
  event_loop([&](){
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    }),
  L_min_max(get_L_min_max(L, x, w_min, w_max)),
  plot(window.geometry(), { .grid = {.color = color::orange} },
        L.to_string(), {}, "w", {}, "L", {}),
  algorithm_event(plot, second_layer),
  algorithm_event_handle(algorithm.event_server().request(algorithm_event, &AlgorithmEvent::callback))
{
  plot.set_xrange({w_min, w_max});
  plot.set_yrange(L_min_max);
  plot.add_to(background_layer, false);

  plot_curve.solve([&L, &x](double w) -> Point { x = w; L.reset_evaluation(); return {w, L.evaluate()}; }, plot.viewport());
  plot.add_bezier_fitter(second_layer, curve_line_style, plot_curve);

  L.reset_evaluation();
}

EnableDrawing::~EnableDrawing()
{
  algorithm_event_handle.cancel();
  window.close();
  event_loop.join();
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  using namespace gradient_descent;

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

    // EXPECTED: learning_rate * derivative was subtracted.
    ASSERT(w == w0 - learning_rate * dLdw);
  }

  //==========================================================================
  Dout(dc::notice, "*** TEST: jump to minimum ***");
  {
    constexpr double w0 = 13.0;

    Algorithm gda(learning_rate, L_max);
    Weight w = w0;
    symbolic::Symbol const& x = symbolic::Symbol::realize("w");
    symbolic::Constant const& c4 = symbolic::Constant::realize(1, 1000000000);
    symbolic::Function const& Lw = symbolic::Function::realize("Lw", (2 - 5 * x + (x^2)) - c4 * (x^4));

//    EnableDrawing enable_drawing(gda, Lw, x, -0.0, 15.0);

    Dout(dc::notice, "Lw = " << fulldef << Lw);

    x = w;
    double dLdw = Lw.derivative(x).evaluate();
    gda(w, Lw.evaluate(), dLdw);                        // 10.3750010985 [one sample, gradient descent]
    ASSERT(utils::almost_equal(static_cast<double>(w), w0 - learning_rate * dLdw, 1e-15));
    Lw.reset_evaluation();
    x = w;
    gda(w, Lw.evaluate(), Lw.derivative(x).evaluate()); // 2.4999957521676026 [jump to vertex]
    ASSERT(utils::almost_equal(static_cast<double>(w), 2.4999957521676026, 1e-15));
    Lw.reset_evaluation();
    x = w;
    // This is the vertex.
    dLdw = Lw.derivative(x).evaluate();
    Dout(dc::notice, "dL/dw = " << dLdw);
    ASSERT(dLdw < 1e-5);
    gda(w, Lw.evaluate(), 0.0);                         // -6231.3483767055095 [best zero]
    // The same value should not be returned twice.
    ASSERT(w != 2.5);
    double best_zero = w;
    ASSERT(best_zero < -6000.0);        // The fourth degree approximation using 13, 10.375 and 2.5 gives something like -6231.
    do
    {
      x = w;
      Lw.reset_evaluation();
      gda(w, Lw.evaluate(), Lw.derivative(x).evaluate()); // -6241.8483809533418 [small step]
                                                          // -13752.619901104652 [continue same direction]
                                                          // 13 [best minimum, opposite direction]
    }
    while (w < best_zero);
    // Since this is the sample that it the furthest away from the vertex and clearly
    // part of the parabolic approximation, that is where we expect to end up again
    // when adding the scale to the (best) minimum.
    ASSERT(w == w0);
  }

  Dout(dc::notice, "Success!");
}
