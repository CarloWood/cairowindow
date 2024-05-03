#include "sys.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
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
#include "debug.h"

class Function
{
 private:
  symbolic::Symbol const& w_ = symbolic::Symbol::realize("w");
  symbolic::Symbol const& a_ = symbolic::Symbol::realize("a");
  symbolic::Symbol const& b_ = symbolic::Symbol::realize("b");
  symbolic::Symbol const& c_ = symbolic::Symbol::realize("c");
  symbolic::Symbol const& half_pi_ = symbolic::Symbol::realize("half_pi");

  symbolic::Expression const& function_ = (a_ + b_ * w_ + c_ * (w_^2)) *
    (half_pi_ * symbolic::atan((a_ + b_ * w_ + c_ * (w_^2) - symbolic::Constant::realize(5.0)) * symbolic::Constant::realize(100.0)) +
     symbolic::Constant::realize(4));
  symbolic::Expression const& derivative_ = function_.derivative(w_);

 public:
  Function()
  {
    a_ = 15.0;
    b_ = 3.1;
    c_ = 0.2;
    half_pi_ = 0.5 * M_PI;
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

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Leaving main()");

  Function L;

  double const w_min = -19.5;
  double const w_max = 4.0;
  int const steps = 100;

  try
  {
    using namespace cairowindow;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Gradient descent of " + L.as_string(), 1200, 900);

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
    }

    // Create and draw plot area.
    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "The function " + L.as_string(), {},
        "w", {},
        "L", {});
    plot.set_xrange({w_min, w_max});
    plot.set_yrange({L_min, L_max});
    plot.add_to(background_layer, false);

    draw::PointStyle point_style({.color_index = 31, .filled_shape = 1});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::LineStyle derivative_line_style({.line_color = color::turquoise, .line_width = 1.0});

    BezierFitter L_fitter;
    L_fitter.solve([&L](double w) -> Point{ return {w, L(w)}; }, {w_min, w_max}, {w_min, L_min, w_max - w_min, L_max - L_min}, 1e-5 * (L_max - L_min));
    auto plot_curve = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(L_fitter));

    auto plot_point = plot.create_point(second_layer, point_style, {50, 0});

    while (true)
    {
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

  Dout(dc::notice, "Entering main()");
}
