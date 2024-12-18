#include "sys.h"
#include "EnableDrawing.h"
#include "Function.h"
#include "cairowindow/EventLoop.h"
#include "gradient_descent2/Algorithm.h"

namespace enable_drawing {
namespace color = cairowindow::color;

EnableDrawing::EnableDrawing(gradient_descent::Algorithm* algorithm, Function const& L,
    double w_min, double w_max, double L_min, double L_max) :
  window("Gradient descent of "
#ifdef CWDEBUG
      + L.to_string(),
#else
      "L",
#endif
      2*600, 2*450),
  background_layer(window.create_background_layer<cairowindow::Layer>(color::silver COMMA_DEBUG_ONLY("background_layer"))),
  second_layer_(window.create_layer<cairowindow::Layer>({} COMMA_DEBUG_ONLY("second_layer"))),
  event_loop([&](){
      // Open window, handle event loop. This must be constructed after the draw stuff, so that it is destructed first!
      // Upon destruction it blocks until the event loop thread finished (aka, the window was closed).
      cairowindow::EventLoop event_loop = window.run();
      event_loop.set_cleanly_terminated();
    }),
  L_min_max(get_L_min_max(L, w_min, w_max)),
  plot_(window.geometry(), { .grid = {.color = color::gray} },
#ifdef CWDEBUG
        L.to_string(),
#else
        "L",
#endif
        {}, "w", {}, "L", {})
        COMMA_DEBUG_ONLY(algorithm_event(plot_, second_layer_),
            algorithm_event_handle(algorithm->event_server().request(algorithm_event, &AlgorithmEvent::callback)))
{
  plot_.set_xrange({w_min, w_max});
  plot_.set_yrange(
      { L_min != std::numeric_limits<double>::max() ? L_min : L_min_max.min(),
        L_max != std::numeric_limits<double>::max() ? L_max : L_min_max.max()});
  plot_.add_to(background_layer, false);

  plot_curve_.solve([&L](double w) -> cairowindow::Point { return {w, L(w)}; }, plot_.viewport());
  plot_.add_bezier_fitter(second_layer_, curve_line_style, plot_curve_);
}

EnableDrawing::~EnableDrawing()
{
#ifdef CWDEBUG
  algorithm_event_handle.cancel();
#endif
  window.close();
  event_loop.join();
}

//static
cairowindow::Range EnableDrawing::get_L_min_max(Function const& L, double w_min, double w_max)
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

} // namespace enable_drawing
