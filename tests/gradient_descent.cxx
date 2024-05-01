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

  int get_zeroes(std::array<double, 2>& zeroes_out) const
  {
    // This can be at most a parabola.
    ASSERT(1 <= coeffients_.size() && coeffients_.size() <= 3);
    if (coeffients_.size() < 3)
    {
      if (coeffients_.size() < 2)
        return 0;
      zeroes_out[0] = -coeffients_[0] / coeffients_[1];
      return 1;
    }

    double D = utils::square(coeffients_[1]) - 4.0 * coeffients_[2] * coeffients_[0];
    if (D < 0.0)
      return 0;
    double delta = std::sqrt(D) / std::abs(2.0 * coeffients_[2]);
    double avg = -coeffients_[1] / (2.0 * coeffients_[2]);
    zeroes_out[0] = avg - delta;
    zeroes_out[1] = avg + delta;
    return utils::almost_equal(zeroes_out[0], zeroes_out[1], 1e-6) ? 1 : 2;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    bool first = true;
    int exponent = 0;
    for (double coeffient : coeffients_)
    {
      if (coeffient != 0.0)
      {
        if (first)
          os << coeffient;
        else if (coeffient > 0.0)
          os << " + " << coeffient;
        else
          os << " - " << -coeffient;
        if (exponent > 0)
        {
          os << ' ' << symbol_name_;
          if (exponent > 1)
            os << '^' << exponent;
        }
        first = false;
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

class Sample
{
 private:
  std::array<double, 4> w_;     // w, w², w³ and w⁴.
  double Lw_;                   // L(w)
  double dLdw_;                 // L'(w)
  cairowindow::plot::Point P_;
  cairowindow::plot::Text P_label_;

 public:
  Sample() = default;
  Sample(double w, double Lw, double dLdw, cairowindow::plot::Point const& point, cairowindow::plot::Text const& label) :
    w_{{ w, w * w, w * w * w, w * w * w *w}}, Lw_(Lw), dLdw_(dLdw), P_(point), P_label_(label) { }

  double w() const { return w_[0]; }
  double pow_w(int exp) const { ASSERT(1 <= exp && exp <= 4); return w_[exp - 1]; }
  double Lw() const { return Lw_; }
  double dLdw() const { return dLdw_; }

  void set_label(cairowindow::plot::Text const& text)
  {
    P_label_ = text;
  }

  cairowindow::Point const& P() const { return P_; }
};

class History
{
 public:
  static constexpr int size = 6;

 private:
  std::array<Sample, size> samples_;
  int prev_ = -1;
  int total_number_of_samples_ = 0;     // The number of samples taken (not necessarily equal to the number that is (still) in the history).
  int current_ = -1;

  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> const& layer_;
  cairowindow::draw::PointStyle const& point_style_;
  cairowindow::draw::TextStyle const& label_style_;

 public:
  History(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer,
      cairowindow::draw::PointStyle const& point_style, cairowindow::draw::TextStyle const& label_style) :
    plot_(plot), layer_(layer), point_style_(point_style), label_style_(label_style) { }

  void add(Sample const& sample)
  {
    current_ = (prev_ + 1) % size;
    samples_[current_] = sample;
    samples_[current_].set_label(plot_.create_text(
          layer_, label_style_({.position = cairowindow::draw::centered_right_of}), sample.P(), std::to_string(total_number_of_samples_)));
    Dout(dc::notice, "Created label " << total_number_of_samples_ << " at w = " << sample.w());
    ++total_number_of_samples_;
  }

  void add(double w, double Lw, double dLdw)
  {
    current_ = (prev_ + 1) % size;
    cairowindow::Point P{w, Lw};
    Sample sample{w, Lw, dLdw,
      plot_.create_point(layer_, point_style_, P),
      plot_.create_text(layer_, label_style_({.position = cairowindow::draw::centered_right_of}), P, std::to_string(total_number_of_samples_))
    };
    Dout(dc::notice, "Appended to history: point " << total_number_of_samples_ << " at w = " << sample.w());
    samples_[current_] = sample;
    ++total_number_of_samples_;
  }

  void append_closest_to(double target)
  {
    ASSERT(total_number_of_samples_ >= 2);
    // The closest sample (in terms of w).
    int closest_index = 0;
    double closest_distance_squared = utils::square(samples_[0].w() - target);
    for (int i = 1; i < std::min(total_number_of_samples_, History::size); ++i)
    {
      double distance_squared = utils::square(samples_[i].w() - target);
      if (distance_squared < closest_distance_squared)
      {
        closest_index = i;
        closest_distance_squared = distance_squared;
      }
    }
    Sample closest = samples_[closest_index];
    advance_prev();
    add(closest);
  }

  static int before(int i) { ASSERT(0 <= i); return (i + size - 1) % size; }

  Sample const& current() const { ASSERT(current_ != -1); ASSERT(current_ < total_number_of_samples_); return samples_[current_]; }
  Sample const& prev() const { ASSERT(prev_ != -1); ASSERT(prev_ < total_number_of_samples_); return samples_[prev_]; }
  Sample const& prev_prev() const { ASSERT(prev_ != -1); ASSERT(prev_ < total_number_of_samples_); return samples_[before(prev_)]; }

  void advance_prev()
  {
    ASSERT(current_ != -1);
    prev_ = current_;
  }

  int total_number_of_samples() const { return total_number_of_samples_; }
};

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Leaving main()");

  Function L;

  double const w_min = 10.0;
  double const w_max = 90.0;
  int const steps = 100;

  double const w_0 = 52.0;
  double learning_rate = 2.0;     // In unit_of(w)^2 / unit_of(L).

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
    L_min = -5.0;
    L_max = 5.0;

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

    // Remember the (most recent) history of samples.
    History history(plot, second_layer, point_style, label_style);

    // Initial value.
    w = w_0;
    double w_delta = 0.0;
    double scale = w_max - w_min;
    Dout(dc::notice, "Initial value: w = " << w);
    int number_of_coef = 0;
    constexpr int down = 1;
    constexpr int up = -1;
    int direction = down;
    std::list<Sample> extremes;
    std::list<Sample>::iterator best_minimum = extremes.end();
    std::list<Sample>::iterator last_extreme = extremes.end();

    while (true)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      double Lw = L(w);
      double dLdw = L.derivative(w);
      history.add(w, Lw, dLdw);
      double prev_w = history.current().w();
      ASSERT(prev_w == w);

      if (number_of_coef < 5)
      {
        ASSERT(history.total_number_of_samples() > 0);
        ASSERT(number_of_coef != 0 || history.total_number_of_samples() == 1);
        ASSERT(number_of_coef != 2 || history.total_number_of_samples() == 2);
        ASSERT(number_of_coef != 3 || history.total_number_of_samples() == 3);
        ASSERT(number_of_coef != 5 || history.total_number_of_samples() >= 4);
        if (history.total_number_of_samples() == 1)
          number_of_coef = 2;
        else if (history.total_number_of_samples() == 2)
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

        Sample const& prev = history.prev();
        double inverse_det = 0.5 / (w - prev.w());
        approximation[1] = inverse_det * 2.0 * (w * prev.dLdw() - prev.w() * dLdw);
        approximation[2] = inverse_det * (dLdw - prev.dLdw());
        parabola_approximation[1] = approximation[1];
        parabola_approximation[2] = approximation[2];
      }
      else
      {
        Sample const& current = history.current();
        Sample const& prev = history.prev();
        Sample const& prev_prev = history.prev_prev();

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
        M <<          1.0,         2.0 *           w,                3.0 *   current.pow_w(2),         4.0 *   current.pow_w(3),
                      1.0,         2.0 *      prev.w(),              3.0 *      prev.pow_w(2),         4.0 *      prev.pow_w(3),
                      1.0,         2.0 * prev_prev.w(),              3.0 * prev_prev.pow_w(2),         4.0 * prev_prev.pow_w(3),
             w - prev.w(), current.pow_w(2) - prev.pow_w(2), current.pow_w(3) - prev.pow_w(3), current.pow_w(4) - prev.pow_w(4);

        Eigen::Vector4d D;
        D << dLdw, prev.dLdw(), prev_prev.dLdw(), Lw - prev.Lw();

        Eigen::Vector4d C = M.colPivHouseholderQr().solve(D);
        approximation[1] = C[0];
        approximation[2] = C[1];
        approximation[3] = C[2];
        approximation[4] = C[3];

        double inverse_det = 0.5 / (w - prev.w());
        parabola_approximation[1] = inverse_det * 2.0 * (w * prev.dLdw() - prev.w() * dLdw);
        parabola_approximation[2] = inverse_det * (dLdw - prev.dLdw());
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

      plot::BezierFitter plot_quotient_curve;

      BezierFitter parabola_approximation_fitter;
      parabola_approximation_fitter.solve([&parabola_approximation](double w) -> Point { return {w, parabola_approximation(w)}; },
          {w_min, w_max}, {w_min, L_min, w_max - w_min, L_max - L_min}, 1e-5 * (L_max - L_min));
      auto plot_parabola_approximation_curve = plot.create_bezier_fitter(second_layer, curve_line_style({.line_color = color::red}),
          std::move(parabola_approximation_fitter));

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // If no extreme exists (or we don't know accurate enough where it is),
      // then `keep_going` will stay true and w is simply adjusted "as usual"
      // by subtracting `learning_rate * dL/dw`.
      bool keep_going = true;
      if (number_of_coef == 3)
      {
        Sample const& prev = history.prev();
        // β = (L'(w₁) - L'(w₀)) / (w₁ - w₀)    [see README.gradient_descent]
        double beta_inverse = (w - prev.w()) / (dLdw - prev.dLdw());
        // If beta is negative, then there is no minimum, only a maximum.
        direction = beta_inverse < 0.0 ? up : down;
        Dout(dc::notice, "direction is set to " << (direction == up ? "up" : "down"));
        ASSERT((approximation[2] < 0.0) == (beta_inverse < 0.0));
        // Set w to the value where the derivative of this parabolic approximation is zero.
        w -= beta_inverse * dLdw;
        keep_going = false;
        // This was the first time we got an idea of the scale at which
        // changes occur. Therefore, use it to set a reasonable learning rate!
        learning_rate = 0.1 * std::abs(beta_inverse);
      }
      else if (number_of_coef == 5)
      {
        Sample const& prev = history.prev();
        double beta_inverse = (w - prev.w()) / (dLdw - prev.dLdw());
        // If beta is negative, then there is no minimum, only a maximum.
        // In that case just keep going in the same direction as before, but accelerating
        // (increase the learning rate with a factor of two).
        if (direction * beta_inverse < 0.0)
          learning_rate *= 2.0;
        else
        {
          // Set w to the value where the derivative of this parabolic approximation is zero.
          double step = -beta_inverse * dLdw;
          w += step;
          keep_going = false;
          // With the new extreme insight, adjust the learning rate to the new scale.
          learning_rate = 0.1 * std::abs(beta_inverse);

          double abs_step = std::abs(step);

          // Are we getting close to the extreme?
          if (abs_step < 0.01 * (w_max - w_min))
          {
            double remainder;
            auto quotient = derivative.long_division(w, remainder);
            Dout(dc::notice, "quotient = " << quotient << " with remainder " << remainder);

            std::array<double, 2> zeroes;
            int number_of_zeroes = quotient.get_zeroes(zeroes);
            if (number_of_zeroes > 1)
              Dout(dc::notice, "with zeroes " << zeroes[0] << " and " << zeroes[1]);
            else if (number_of_zeroes == 1)
              Dout(dc::notice, "with one zero at " << zeroes[0]);
            else
              Dout(dc::notice, "with no zeroes!");

            plot_quotient_curve.solve([&quotient](double w) -> Point { return {w, 10.0 * quotient(w)}; },
                {w_min, w_max}, {w_min, L_min, w_max - w_min, L_max - L_min}, 1e-5 * (L_max - L_min));
            plot.add_bezier_fitter(second_layer, curve_line_style({.line_color = color::blue}), plot_quotient_curve);

            if (number_of_zeroes == 2)
            {
              scale = std::min(scale, zeroes[1] - zeroes[0]);
              Dout(dc::notice, "scale was set to " << scale);
            }

            // Did we reach the (local) extreme?
            if (number_of_zeroes == 2 && abs_step < 0.01 * scale)
            {
              Dout(dc::notice, (direction == up ? "Maximum" : "Minimum") << " reached!");

              // Reset scale.
              scale = w_max - w_min;
              Dout(dc::notice, "scale reset to " << scale);

              // Change direction.
              direction = -direction;
              Dout(dc::notice, "direction is set to " << (direction == up ? "up" : "down") << "; w_delta = " << w_delta);

              int best_zero = (number_of_zeroes == 2 &&
                  (w_delta > 0.0 || (w_delta == 0.0 && approximation(zeroes[1]) < approximation(zeroes[0])))) ? 1 : 0;
              w = zeroes[best_zero];

              // Now that the decision on which direction (left/right) we explore is taken,
              // store that decision (as the sign of w_delta);
              if (w_delta == 0.0)
              {
                w_delta = w - history.current().w();
                Dout(dc::notice, "w_delta --> " << w_delta);
              }

              // Store this extreme.
              std::list<Sample>::iterator new_extreme;
              if (w_delta > 0.0)
                new_extreme = extremes.insert(extremes.end(), history.current());
              else
                new_extreme = extremes.insert(extremes.begin(), history.current());

              if (direction == up)
              {
                if (best_minimum == extremes.end() || best_minimum->Lw() > new_extreme->Lw())
                {
                  best_minimum = new_extreme;
                  Dout(dc::notice, "best_minimum set to " << best_minimum->w() << " / " << best_minimum->Lw());
                }
                if (new_extreme != best_minimum)
                {
                  // The new minimum isn't better than what we found already. Stop going into this direction.
                  break;
                }
              }

              // Re-add old sample again to the history; add one closest to the target.
              history.append_closest_to(w);
              history.advance_prev();
              // Add the target.
              Lw = L(w);
              dLdw = L.derivative(w);
              history.add(w, Lw, dLdw);
              // Add one more, using the learning rate.
              keep_going = true;
            }
          }
        }
      }

      if (keep_going)
      {
        // There is no new extreme insight yet; just keep going (up or) down hill.
        if (direction == down)
          w -= learning_rate * dLdw;
        else
          w += learning_rate * dLdw;
      }

      Dout(dc::notice, history.total_number_of_samples() << ": " << prev_w << " --> " << w <<
          " (learning rate is now " << learning_rate << ")");
      history.advance_prev();

      // Block until a redraw is necessary (for example because the user moved a draggable object,
      // or wants to print the current drawing) then go to the top of loop for a redraw.
      if (!window.handle_input_events())
        break;          // Program must be terminated.

//      if (utils::almost_equal(history.prev().w(), w, 1e-6))
//        break;
    }

    Dout(dc::notice, "Found global minimum " << best_minimum->Lw() << " at w = " << best_minimum->w());
    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Entering main()");
}
