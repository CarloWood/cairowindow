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
#include "debug.h"

class Function
{
 private:
  symbolic::Symbol const& w_ = symbolic::Symbol::realize("w");
  symbolic::Symbol const& a_ = symbolic::Symbol::realize("a");
  symbolic::Symbol const& b_ = symbolic::Symbol::realize("b");
  symbolic::Symbol const& c_ = symbolic::Symbol::realize("c");
  symbolic::Symbol const& d_ = symbolic::Symbol::realize("d");
  symbolic::Symbol const& e_ = symbolic::Symbol::realize("e");

  //symbolic::Expression const& function_ = a_ + b_ * w_ + c_ * (w_^2) + d_ * (w_^3) + e_ * (w_^4);
  symbolic::Expression const& function_ = symbolic::sin(a_ * w_);
  symbolic::Expression const& derivative_ = function_.derivative(w_);

 public:
  Function()
  {
    //    4        3
    //  w    130⋅w          2
    //  ── - ────── + 2200⋅w  - 32000⋅w + 10000
    //  4      3

    a_ = 0.1; //10000.0;
    b_ = -32000.0;
    c_ = 2200.0;
    d_ = -130.0 / 3.0;
    e_ = 0.25;
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

class Data
{
 private:
  std::array<double, 4> w_;     // w, w², w³ and w⁴.
  double Lw_;                   // L(w)
  double dLdw_;                 // L'(w)

 public:
  Data() = default;
  Data(double w, double Lw, double dLdw) : w_{{ w, w * w, w * w * w, w * w * w *w}}, Lw_(Lw), dLdw_(dLdw) { }

  double w() const { return w_[0]; }
  double pow_w(int exp) const { ASSERT(1 <= exp && exp <= 4); return w_[exp - 1]; }
  double Lw() const { return Lw_; }
  double dLdw() const { return dLdw_; }
};

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Leaving main()");

  Function L;

  double const w_min = -20.0;
  double const w_max = 100.0;
  int const steps = 100;

  double const w_0 = 20.0;
  double learning_rate = 0.1;     // In unit_of(w)^2 / unit_of(L).

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

    std::vector<Point> curve_points;
    double w = w_min;
    double delta_w = (w_max - w_min) / (steps - 1);
    double L_min = L(w_min);
    double L_max = L_min;
    for (int i = 0; i < steps; ++i)
    {
      double val = L(w);
      curve_points.emplace_back(w, val);
      L_min= std::min(L_min, val);
      L_max= std::max(L_max, val);
      w += delta_w;
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

    // Remember the prev three measurements.
    std::array<Data, 3> history;
    int prev = -1;

    // Initial value.
    w = w_0;
    int number_of_coef;

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      int cur = (prev + 1) % history.size();
      history[cur] = Data{w, L(w), L.derivative(w)};
      if (number_of_coef < 5)
      {
        if (cur == 0)
          number_of_coef = 2;
        else if (cur == 1)
          number_of_coef = 3;
        else
          number_of_coef = 5;
      }

      // Convenience variables.
      double const Lw = history[cur].Lw();
      double const dLdw = history[cur].dLdw();

      // Create a point P.
      auto plot_P = plot.create_point(second_layer, point_style, {w, Lw});
      auto P_label = plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), plot_P, "P");

      std::vector<double> coef(number_of_coef);
      if (number_of_coef == 2)
      {
        // If we have just one point, then the approximation is a linear function:
        //
        // A(w) = coef[0] + L'(w) w
        coef[1] = dLdw;
      }
      else if (number_of_coef == 3)
      {
        // If we have two points, the approximation is a parabola:
        //
        //   A(w) = coef[0] + b w + c w²
        //
        // for which we determined the value of the derivative at two points:
        //
        //   L'(w₀) = b + 2c w₀
        //   L'(w₁) = b + 2c w₁
        //
        // In matrix form:
        //
        //   ⎡1 2w₀⎤ ⎡b⎤   ⎡L'(w₀)⎤
        //   ⎣1 2w₁⎦ ⎣c⎦ = ⎣L'(w₁)⎦
        //
        // from which follows
        //
        // ⎡b⎤        1     ⎡2w₁ -2w₀⎤⎡L'(w₀)⎤        1      ⎡2w₁ L'(w₀) - 2w₀ L'(w₁)⎤
        // ⎣c⎦ = ---------- ⎣-1   1  ⎦⎣L'(w₁)⎦ = ----------- ⎣    L'(w₁) -     L'(w₀)⎦
        //        2w₁ - 2w₀                      2 (w₁ - w₀)

        double inverse_det = 0.5 / (w - history[prev].w());
        coef[1] = inverse_det * 2.0 * (w * history[prev].dLdw() - history[prev].w() * dLdw);
        coef[2] = inverse_det * (dLdw - history[prev].dLdw());
      }
      else
      {
        // If we have (at least) three points, the approximation is a fourth degree parabola:
        //
        //   A(w) = a + b w + c w² + d w³ + e w⁴
        //
        // for which we determined the value of the derivative at three points, L'(w₀), L'(w₁) and L'(w₂).
        //
        // The matrix form for the coefficients of the polynomial then becomes:
        //
        //   ⎡  1      2w₂      3w₂²     4w₂³ ⎤ ⎡b⎤   ⎡    L'(w₂)   ⎤
        //   ⎢  1      2w₁      3w₁²     4w₁³ ⎥ ⎢c⎥   ⎢    L'(w₁)   ⎥
        //   ⎢  1      2w₀      3w₀²     4w₀³ ⎥ ⎢d⎥ = ⎢    L'(w₀)   ⎥
        //   ⎣w₂-w₁  w₂²-w₁²  w₂³-w₁³  w₂⁴-w₁⁴⎦ ⎣e⎦   ⎣L(w₂) - L(w₁)⎦
        //
        Eigen::Matrix4d M;
        int prev_prev = (prev + (history.size() - 1)) % history.size();

        M << 1.0, 2.0 * w, 3.0 * history[cur].pow_w(2), 4.0 * history[cur].pow_w(3),
             1.0, 2.0 * history[prev].w(), 3.0 * history[prev].pow_w(2), 4.0 * history[prev].pow_w(3),
             1.0, 2.0 * history[prev_prev].w(), 3.0 * history[prev_prev].pow_w(2), 4.0 * history[prev_prev].pow_w(3),
             w - history[prev].w(), history[cur].pow_w(2) - history[prev].pow_w(2), history[cur].pow_w(3) - history[prev].pow_w(3), history[cur].pow_w(4) - history[prev].pow_w(4);

        Eigen::Vector4d D;
        D << dLdw, history[prev].dLdw(), history[prev_prev].dLdw(), Lw - history[prev].Lw();

        Eigen::Vector4d C = M.colPivHouseholderQr().solve(D);
        coef[1] = C[0];
        coef[2] = C[1];
        coef[3] = C[2];
        coef[4] = C[3];
      }
      auto approximation = [&](double w) -> Point
      {
        double result = 0.0;
        for (double c : std::ranges::reverse_view(coef))
          result = w * result + c;
        return {w, result};
      };
      coef[0] = Lw - approximation(w).y();

      BezierFitter approximation_fitter;
      approximation_fitter.solve(approximation, {w_min, w_max}, {w_min, L_min, w_max - w_min, L_max - L_min}, 1e-5 * (L_max - L_min));
      auto plot_approximation_curve = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::teal}),
          std::move(approximation_fitter));

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      if (number_of_coef == 3)
      {
        double new_learning_rate = (w - history[prev].w()) / (dLdw - history[prev].dLdw());
        // If the new learning rate is negative, then there is no minimum, only a maximum.
        // In that case just keep going in the same direction as before, but accelerating
        // (increase the learning rate with a factor of two).
        if (new_learning_rate < 0.0)
          learning_rate *= 2.0;
        else
          learning_rate = new_learning_rate;
      }
      else if (number_of_coef == 5)
      {
        //FIXME: implement this.
        double new_learning_rate = (w - history[prev].w()) / (dLdw - history[prev].dLdw());
        // If the new learning rate is negative, then there is no minimum, only a maximum.
        // In that case just keep going in the same direction as before, but accelerating
        // (increase the learning rate with a factor of two).
        if (new_learning_rate < 0.0)
          learning_rate *= 2.0;
        else
          learning_rate = new_learning_rate;
      }

      // Adjust w.
      w = w - learning_rate * dLdw;

      prev = cur;
      Dout(dc::notice, history[prev].w() << " --> " << w << " (learning rate is now " << learning_rate << ")");

      // Block until a redraw is necessary (for example because the user moved a draggable object,
      // or wants to print the current drawing) then go to the top of loop for a redraw.
      if (!window.handle_input_events())
        break;          // Program must be terminated.

      if (utils::almost_equal(history[prev].w(), w, 1e-6))
        break;
    }

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Entering main()");
}
