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

#if 0
  Expression const& polynomial = a + b * w + c * (w^2) + d * (w^3) + e * (w^4) + f * (w^5);
  Dout(dc::notice, "polynomial = " << polynomial);
  a = -2.0;
  b = 3.0;
  c = -0.2;
  d = 0.4;
  e = -0.01;
  f = 0.002;
#endif

using utils::has_print_on::operator<<;
class Polynomial
{
 private:
  std::vector<double> coeffients_;
#ifdef CWDEBUG
  std::string symbol_name_;
#endif

 public:
  // Create a polynomial of at most degree `number_of_coeffients - 1`, starting
  // with all coeffients set to zero.
  Polynomial(int number_of_coeffients COMMA_CWDEBUG_ONLY(std::string const& symbol_name)) :
    coeffients_(number_of_coeffients) COMMA_CWDEBUG_ONLY(symbol_name_(symbol_name))
  {
    ASSERT(coeffients_.size() == number_of_coeffients);
    for (int i = 0; i < number_of_coeffients; ++i)
      ASSERT(coeffients_[i] == 0.0);
  }

  double operator[](int i) const { ASSERT(0 <= i && i < coeffients_.size()); return coeffients_[i]; }
  double& operator[](int i) { ASSERT(0 <= i && i < coeffients_.size()); return coeffients_[i]; }

  double operator()(double w)
  {
    double result = 0.0;
    for (double coeffient : std::ranges::reverse_view(coeffients_))
      result = w * result + coeffient;
    return result;
  };

  Polynomial derivative() const
  {
    Polynomial result(coeffients_.size() - 1 COMMA_CWDEBUG_ONLY(symbol_name_));
    for (int i = 1; i < coeffients_.size(); ++i)
      result[i - 1] = coeffients_[i] * i;
    return result;
  }

  // Return the division of this Polynomial by the factor (w - z).
  Polynomial long_division(double z, double& remainder) const
  {
    // f(w) = 3 * w^3 +      5 * w^2 - 4 * w + 10.
    //        3 * w^3 + (-2)*3 * w^2
    //      - ------------------------------------
    //                      11 * w^2 -     4 * w + 10.
    //                      11 * w^2 + (-2)*11 w
    //                    - --------------------------
    //                                    18 * w +     10.
    //                                    18 * w + (-2)18
    //                                  - ---------------
    //                                                 46
    // Divide by (w - 2)
    // 3 * w^2 + 11 * w + 18

    // NOTICE        : 10 + -4 w + 5 w^2 + 3 w^3
    // (w - 2)(3 w^2 + 11 w + 18) = 3 w^3 + 5 w^2 - 4 w - 36

    if (coeffients_.size() < 2)
    {
      ASSERT(coeffients_.size() == 1);
      remainder = coeffients_[0];
      return {1 COMMA_CWDEBUG_ONLY(symbol_name_)};
    }
    Polynomial result(coeffients_.size() - 1 COMMA_CWDEBUG_ONLY(symbol_name_));
    result[coeffients_.size() - 2] = coeffients_[coeffients_.size() - 1];
    for (int i  = coeffients_.size() - 2; i > 0; --i)
      result[i - 1] = coeffients_[i] + z * result[i];
    remainder = coeffients_[0] + z * result[0];
    return result;
  }

  std::pair<double, double> zeroes() const
  {
    // This can be at most a parabola.
    ASSERT(1 <= coeffients_.size() && coeffients_.size() <= 3);
    if (coeffients_.size() < 3)
    {
      if (coeffients_.size() < 2)
        return {0.0, 0.0};
      double zero = -coeffients_[0] / coeffients_[1];
      return {zero, zero};
    }

    double D = utils::square(coeffients_[1]) - 4.0 * coeffients_[2] * coeffients_[0];
    ASSERT(D >= 0.0);
    double delta = std::sqrt(D) / std::abs(2.0 * coeffients_[2]);
    double avg = -coeffients_[1] / (2.0 * coeffients_[2]);
    return {avg - delta, avg + delta};
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    char const* sep = "";
    int exponent = 0;
    for (double coeffient : coeffients_)
    {
      if (coeffient != 0.0)
      {
        os << sep << coeffient;
        if (exponent > 0)
        {
          os << ' ' << symbol_name_;
          if (exponent > 1)
            os << '^' << exponent;
        }
        sep = " + ";
      }
      ++exponent;
    }
  }
#endif
};

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
    b_ = 0.15;
    c_ = 57.0;
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

class Data
{
 private:
  std::array<double, 4> w_;     // w, w², w³ and w⁴.
  double Lw_;                   // L(w)
  double dLdw_;                 // L'(w)
  cairowindow::plot::Point P_;
  cairowindow::plot::Text P_label_;

 public:
  Data() = default;
  Data(double w, double Lw, double dLdw, cairowindow::plot::Point const& point, cairowindow::plot::Text const& label) :
    w_{{ w, w * w, w * w * w, w * w * w *w}}, Lw_(Lw), dLdw_(dLdw), P_(point), P_label_(label) { }

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

  double const w_0 = 32.0;
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
    int number_of_coef = 0;
    int n = 0;

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      int cur = (prev + 1) % history.size();
      double const Lw = L(w);
      double const dLdw = L.derivative(w);
      Point P{w, Lw};
      history[cur] = Data{w, Lw, dLdw,
        plot.create_point(second_layer, point_style, P),
        plot.create_text(second_layer, label_style({.position = draw::centered_right_of}), P, std::to_string(n))
      };
      ++n;
      if (number_of_coef < 5)
      {
        if (cur == 0)
          number_of_coef = 2;
        else if (cur == 1)
          number_of_coef = 3;
        else
          number_of_coef = 5;
      }

      Polynomial parabola_approximation(3 COMMA_CWDEBUG_ONLY("w"));

      Polynomial approximation(number_of_coef COMMA_CWDEBUG_ONLY("w"));
      if (number_of_coef == 2)
      {
        // If we have just one point, then the approximation is a linear function:
        //
        // A(w) = coef[0] + L'(w) w
        approximation[1] = dLdw;
        parabola_approximation[1] = dLdw;
      }
      else if (number_of_coef == 3)
      {
        // If we have two points, the approximation is a parabola:
        //
        //   A(w) = approximation[0] + b w + c w²
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
        approximation[1] = inverse_det * 2.0 * (w * history[prev].dLdw() - history[prev].w() * dLdw);
        approximation[2] = inverse_det * (dLdw - history[prev].dLdw());
        parabola_approximation[1] = approximation[1];
        parabola_approximation[2] = approximation[2];
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
        approximation[1] = C[0];
        approximation[2] = C[1];
        approximation[3] = C[2];
        approximation[4] = C[3];

        double inverse_det = 0.5 / (w - history[prev].w());
        parabola_approximation[1] = inverse_det * 2.0 * (w * history[prev].dLdw() - history[prev].w() * dLdw);
        parabola_approximation[2] = inverse_det * (dLdw - history[prev].dLdw());
      }
      approximation[0] = Lw - approximation(w);
      parabola_approximation[0] = Lw - parabola_approximation(w);
      Dout(dc::notice, "approximation = " << approximation);
      auto derivative = approximation.derivative();
      Dout(dc::notice, "derivative = " << derivative);

      BezierFitter approximation_fitter;
      approximation_fitter.solve([&approximation](double w) -> Point { return {w, approximation(w)}; },
          {w_min, w_max}, {w_min, L_min, w_max - w_min, L_max - L_min}, 1e-5 * (L_max - L_min));
      auto plot_approximation_curve = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::teal}),
          std::move(approximation_fitter));

      BezierFitter derivative_fitter;
      derivative_fitter.solve([&derivative](double w) -> Point { return {w, derivative(w)}; },
          {w_min, w_max}, {w_min, L_min, w_max - w_min, L_max - L_min}, 1e-5 * (L_max - L_min));
      auto plot_derivative_curve = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::magenta}),
          std::move(derivative_fitter));

      double remainder;
      auto quotient = derivative.long_division(w, remainder);
      Dout(dc::notice, "quotient = " << quotient);

      auto zeroes = quotient.zeroes();
      Dout(dc::notice, "with zeroes " << zeroes.first << " and " << zeroes.second);

      BezierFitter quotient_fitter;
      quotient_fitter.solve([&quotient](double w) -> Point { return {w, 10.0 * quotient(w)}; },
          {w_min, w_max}, {w_min, L_min, w_max - w_min, L_max - L_min}, 1e-5 * (L_max - L_min));
      auto plot_quotient_curve = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::blue}),
          std::move(quotient_fitter));

      BezierFitter parabola_approximation_fitter;
      parabola_approximation_fitter.solve([&parabola_approximation](double w) -> Point { return {w, parabola_approximation(w)}; },
          {w_min, w_max}, {w_min, L_min, w_max - w_min, L_max - L_min}, 1e-5 * (L_max - L_min));
      auto plot_parabola_approximation_curve = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::red}),
          std::move(parabola_approximation_fitter));

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      bool keep_going = true;
      if (number_of_coef == 3)
      {
        // β = (L'(w₁) - L'(w₀)) / (w₁ - w₀)    [see README.gradient_descent]
        double beta_inverse = (w - history[prev].w()) / (dLdw - history[prev].dLdw());
        // If beta is negative, then there is no minimum, only a maximum.
        // In that case just keep going in the same direction as before, but accelerating
        // (increase the learning rate with a factor of two).
        if (beta_inverse < 0.0)
          learning_rate *= 2.0;
        else
        {
          // We see a minimum. Hit the brakes!
          learning_rate = beta_inverse;
          // Set w to the value where the derivative of this parabolic approximation is zero.
          w -= beta_inverse * dLdw;
          keep_going = false;
        }
      }
      else if (number_of_coef == 5)
      {
        //FIXME: implement this.
        double beta_inverse = (w - history[prev].w()) / (dLdw - history[prev].dLdw());
        // If beta is negative, then there is no minimum, only a maximum.
        // In that case just keep going in the same direction as before, but accelerating
        // (increase the learning rate with a factor of two).
        if (beta_inverse < 0.0)
          learning_rate *= 2.0;
        else
        {
          // We see a minimum. Hit the brakes!
          learning_rate = beta_inverse;
          // Set w to the value where the derivative of this parabolic approximation is zero.
          w -= beta_inverse * dLdw;
          keep_going = false;

          // Did we reach a (local) minimum?
          if (std::abs(w - history[prev].w()) < 0.01 * (zeroes.second - zeroes.first))
          {
            Dout(dc::notice, "Minimum reached!");
          }
        }
      }

      if (keep_going)
      {
        // There is no new minimum known; just keep going down hill.
        w = w - learning_rate * dLdw;
      }

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
