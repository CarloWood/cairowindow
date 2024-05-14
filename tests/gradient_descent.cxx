#include "sys.h"
#include "Polynomial.h"
#include "QuadraticPolynomial.h"
#include "Sample.h"
#include "Scale.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/symbolic/symbolic.h"
#include "utils/square.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
#include "utils/AIRefCount.h"
#include <Eigen/Dense>
#include <sstream>
#include <thread>
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

#if 0
class Function
{
 public
  static constexpr double w_0 = 12.0;
  static constexpr double w_min = -20.0;
  static constexpr double w_max = 80.0;

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
#else
constexpr double a2 = 14.6;
constexpr double b2 = 3.15;
constexpr double c2 = 0.451;

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
    a2_ = a2;
    b2_ = b2;
    c2_ = c2;
  }

  void set_a1(double a1)
  {
    a1_ = a1;
  }

  void set_b1(double b1)
  {
    b1_ = b1;
  }

  void set_c1(double c1)
  {
    c1_ = c1;
  }

  void set_a2(double a2)
  {
    a2_ = a2;
  }

  void set_b2(double b2)
  {
    b2_ = b2;
  }

  void set_c2(double c2)
  {
    c2_ = c2;
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

// A Sample including plot objects: a point and a label.
class Sample : public gradient_descent::Sample
{
 private:
  mutable cairowindow::plot::Point P_;
  mutable cairowindow::plot::Text P_label_;

 public:
  // Create an uninitialized Sample, required for the History.
  Sample() = default;

  // Constructor.
  Sample(Weight w, double Lw, double dLdw,
      cairowindow::plot::Point const& point, cairowindow::plot::Text const& label) :
    gradient_descent::Sample(w, Lw, dLdw), P_(point), P_label_(label) { }

  void set_point(cairowindow::plot::Point const& P) const
  {
    P_ = P;
  }

  void set_label(cairowindow::plot::Text const& text) const
  {
    P_label_ = text;
  }

  cairowindow::Point const& P() const { return P_; }
  cairowindow::Text const& label() const { return P_label_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << label().text() << " (at " << w_ << ")";
  }
#endif
};

class Scale : public gradient_descent::Scale
{
 public:
  static cairowindow::draw::ConnectorStyle const s_indicator_style;

 private:
  // Used to visualize the Scale:
  cairowindow::plot::Plot* plot_ = nullptr;
  boost::intrusive_ptr<cairowindow::Layer> layer_;
  cairowindow::plot::Connector plot_indicator_;
  cairowindow::plot::Text plot_scale_text_;
  cairowindow::plot::Line plot_vertical_line_through_w_;
  cairowindow::plot::Line plot_vertical_line_through_v_;

  // Temporary curves (used while developing this class).
  cairowindow::plot::BezierFitter plot_old_parabola_;

 private:
  friend class gradient_descent::Scale;
  Scale(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer) : plot_(&plot), layer_(layer) { }

 public:
  void draw_indicators()
  {
    ASSERT(has_sample_);

    double x1 = parabola_.vertex_x();
    double x2 = edge_sample_w_;
    double scale_y = plot_->yrange().min() + 0.5 * plot_->yrange().size();
    plot_indicator_ = cairowindow::plot::Connector{{x1, scale_y}, {x2, scale_y},
        cairowindow::Connector::open_arrow, cairowindow::Connector::open_arrow};
    plot_->add_connector(layer_, s_indicator_style, plot_indicator_);
    plot_scale_text_ = plot_->create_text(layer_, {{.position = cairowindow::draw::centered_above, .offset = 2.0}},
        cairowindow::Point{(x1 + x2) / 2, scale_y}, "scale");

    plot_vertical_line_through_v_ = cairowindow::plot::Line{{parabola_.vertex_x(), scale_y}, cairowindow::Direction::up};
    plot_->add_line(layer_, s_indicator_style({.dashes = {3.0, 3.0}}), plot_vertical_line_through_v_);

    plot_vertical_line_through_w_ = cairowindow::plot::Line{{edge_sample_w_, scale_y}, cairowindow::Direction::up};
    plot_->add_line(layer_, s_indicator_style, plot_vertical_line_through_w_);
  }

  void erase_indicators()
  {
    // Erase old drawings.
    plot_indicator_.reset();
    plot_scale_text_.reset();
    plot_vertical_line_through_w_.reset();
    plot_vertical_line_through_v_.reset();
    plot_old_parabola_.reset();
  }

  // For use with draw_old_parabola.
  math::QuadraticPolynomial get_parabola() const
  {
    return parabola_;
  }

  // Draws and caches a plot of the old parabola; which has to be passed as an argument.
  void draw_old_parabola(math::QuadraticPolynomial const& old_parabola)
  {
    // Draw the old parabola, for debugging purposes.
    using namespace cairowindow;
    plot_old_parabola_.solve([&](double w) -> Point { return {w, old_parabola(w)}; }, plot_->viewport());
    plot_->add_bezier_fitter(layer_, {{.line_color = color::light_red, .line_width = 1.0}}, plot_old_parabola_);
  }
};

class SampleMinimum : public Sample
{
 private:
  boost::intrusive_ptr<Scale> scale_;
  double energy_;
  int done_;

 public:
  SampleMinimum(Sample const& minimum, boost::intrusive_ptr<Scale> const& scale, double energy, int hdirection) :
    Sample(minimum), scale_(scale), energy_(energy), done_(hdirection == -1 ? 1 : 2) { }

  boost::intrusive_ptr<Scale> scale(int hdirection)
  {
    int done_flag = hdirection == -1 ? 1 : 2;
    ASSERT((done_ & done_flag) == 0);
    done_ |= done_flag;
    return scale_;
  }
  double energy() const { return energy_; }

  bool done() const { return done_ == 3; }
};

namespace gradient_descent {

template<ConceptSample T>
class History
{
 public:
  static constexpr int size = 9;

 private:
  std::array<T, size> samples_;
  int current_ = -1;                    // Index of the last Sample that was added.
  int prev_ = -1;                       // Index to the Sample that was added before current.
  int old_samples_ = 0;                 // Samples elsewhere that should not be used for fitting a polynomial.
  int total_number_of_samples_ = 0;     // The number of samples taken (not necessarily equal to the number that is (still) in the history).

 public:
  T const& add(double w, double Lw, double dLdw, boost::intrusive_ptr<Scale> const& scale)
  {
    // If the new sample is very close to the current one, don't add it, just replace the current sample.
    bool almost_equal = (current_ != -1 && std::abs(samples_[current_].w() - w) < 0.001 * scale->or_zero());
    if (almost_equal)
    {
      Dout(dc::notice, "Replacing history point " << (total_number_of_samples_ - 1) << " at " << samples_[current_].w() << " --> " << w);
    }
    else
    {
      Dout(dc::notice, "Appending to history: point " << total_number_of_samples_ << " at w = " << w);

      prev_ = current_;
      current_ = (prev_ + 1) % size;
      ++total_number_of_samples_;
    }

    samples_[current_].set_values(w, Lw, dLdw);
    return samples_[current_];
  }

  void reset()
  {
    old_samples_ = total_number_of_samples_;
  }

  int relevant_samples() const
  {
    return total_number_of_samples_ - old_samples_;
  }

  static int before(int i) { ASSERT(0 <= i); return (i + size - 1) % size; }

  T const& current() const { ASSERT(current_ != -1); ASSERT(current_ < total_number_of_samples_); return samples_[current_]; }
  T const& prev() const { ASSERT(prev_ != -1); ASSERT(prev_ < total_number_of_samples_); return samples_[prev_]; }
  T const& prev_prev() const { ASSERT(prev_ != -1); ASSERT(prev_ < total_number_of_samples_); return samples_[before(prev_)]; }

  int total_number_of_samples() const { return total_number_of_samples_; }
};

} // namespace gradient_descent

template<gradient_descent::ConceptSample T>
class History : public gradient_descent::History<T>
{
 private:
  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> const& layer_;
  cairowindow::draw::PointStyle const& point_style_;
  cairowindow::draw::TextStyle const& label_style_;

 public:
  History(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer,
      cairowindow::draw::PointStyle const& point_style, cairowindow::draw::TextStyle const& label_style) :
    plot_(plot), layer_(layer), point_style_(point_style), label_style_(label_style) { }
};

//static
cairowindow::draw::ConnectorStyle const Scale::s_indicator_style{{.line_width = 1}};

// A class describing the amount of "kinetic energy" (elsewhere in literature more often known as "momentum")
// that we have, which determines the maximum overshoot uphill that we allow. Having a certain momentum or
// kinetic energy helps in overcoming small bumbs and not get stuck in the first local minimum that we encounter.
class KineticEnergy
{
  static constexpr double friction = 0.0513;    // (1 - friction)^2 = 0.9
  static cairowindow::draw::ConnectorStyle const s_indicator_style;

 private:
  double max_Lw_;       // The maximum height that could be reached with the current amount of kinetic energy.
  double Lw_;           // The current height. The kinetic energy is proportional to max_Lw_ - Lw_.

  // Used to visualize the KineticEnergy:
  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> layer_;
  cairowindow::plot::Line plot_horizontal_line_;
  cairowindow::plot::Text plot_energy_text_;

 private:
  void draw_indicators()
  {
    plot_horizontal_line_ = cairowindow::plot::Line{{0.0, max_Lw_}, cairowindow::Direction::right};
    plot_.add_line(layer_, s_indicator_style, plot_horizontal_line_);

    plot_energy_text_ = plot_.create_text(layer_, {{.position = cairowindow::draw::centered_above, .offset = 2.0}},
        cairowindow::Point{0.5 * (plot_.xrange().min() + plot_.xrange().max()), max_Lw_}, "energy");
  }

  double max_Lw(double new_Lw) const
  {
    // Reduce the maximum height with a fraction of the (vertical) distance traveled.
    return max_Lw_ - friction * std::abs(new_Lw - Lw_);
  }

 public:
  KineticEnergy(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer, double Lw) :
    plot_(plot), layer_(layer), max_Lw_(Lw), Lw_(Lw) { }

  void set(double max_Lw, double Lw)
  {
    DoutEntering(dc::notice, "KineticEnergy::set(" << max_Lw << ", " << Lw << ")");

    max_Lw_ = max_Lw;
    Lw_ = Lw;

    draw_indicators();
    Dout(dc::notice, "max_Lw_ set to " << max_Lw_ << "; kinetic energy is now " << (max_Lw_ - Lw_));
  }

  double energy() const
  {
    return max_Lw_;             // This is the total energy.
  }

  // Returns true upon success. If false is returned the update is rejected and should not take place!
  bool maybe_update(double new_Lw)
  {
    DoutEntering(dc::notice, "KineticEnergy::maybe_update(" << new_Lw << ")");

    // Reduce the maximum height with a fraction of this distance.
    double new_max_Lw = max_Lw(new_Lw);

    // Check if this is a request to go higher than the maximum height.
    if (new_Lw > new_max_Lw)
    {
      Dout(dc::notice, "Rejected because max_Lw_ = " << max_Lw_ << ", Lw_ = " << Lw_ << "; new_Lw = " << new_Lw << " > " <<
          "max_Lw_ - friction * std::abs(new_Lw - Lw_) = " << new_max_Lw);
      return false;
    }

    max_Lw_ = new_max_Lw;
    Lw_ = new_Lw;

    draw_indicators();
    Dout(dc::notice, "max_Lw_ set to " << max_Lw_ << "; kinetic energy is now " << (max_Lw_ - Lw_));
    return true;
  }

  void update(double new_Lw)
  {
    DoutEntering(dc::notice, "KineticEnergy::update(" << new_Lw << ")");

    // Reduce total energy (the maximum height) as usual, but if we need
    // more energy than we have, just reset the total energy to zero.
    max_Lw_ = std::max(new_Lw, max_Lw(new_Lw));
    Lw_ = new_Lw;

    draw_indicators();
    Dout(dc::notice, "max_Lw_ set to " << max_Lw_ << "; kinetic energy is now " << (max_Lw_ - Lw_));
  }

#if CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{max height:" << max_Lw_ << ", energy:" << (max_Lw_ - Lw_) << "}";
  }
#endif
};

//static
cairowindow::draw::ConnectorStyle const KineticEnergy::s_indicator_style{{.line_width = 1}};

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  Function L;

  double const w_0 = Function::w_0;
  double const w_min = Function::w_min;
  double const w_max = Function::w_max;

  int const steps = 100;
  double const learning_rate = 0.1;     // In unit_of(w)^2 / unit_of(L).

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
//      L_min = -1.25;
//      L_max = 60.0;
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

    BezierFitter L_fitter([&L](double w) -> Point { return {w, L(w)}; }, plot.viewport());
    auto plot_curve = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(L_fitter));

    // Remember the (most recent) history of samples.
    History<Sample> history(plot, second_layer, point_style, label_style);

    // Initial value.
#if 0
    double Lw_0 = L(w_0);
    Sample new_sample(w_0, Lw_0, L.derivative(w_0),
      plot.create_point(second_layer, point_style, cairowindow::Point{w_0, Lw_0}),
      plot.create_text(second_layer, label_style({.position = cairowindow::draw::centered_below}), cairowindow::Point{w_0, Lw_0}, "0"));
#endif
    gradient_descent::Weight w(w_0);
    Dout(dc::notice, "Initial value of new_sample; w = " << w);
    double current_learning_rate = learning_rate;
    auto scale = Scale::create<Scale>(plot, second_layer);
    KineticEnergy energy(plot, second_layer, L_max);
    Dout(dc::notice, "Initial height is " << energy);
    constexpr int unknown = 0;
    constexpr int down = 1;
    constexpr int up = -1;
    constexpr int left = -1;
    constexpr int right = 1;
    int vdirection = down;
    int hdirection = unknown;
    std::list<SampleMinimum> extremes;
    std::list<SampleMinimum>::iterator best_minimum = extremes.end();
    std::list<SampleMinimum>::iterator last_extreme = extremes.end();

    plot::BezierFitter plot_approximation_curve;
    plot::BezierFitter plot_derivative_curve;
    plot::BezierFitter plot_quotient_curve;
    plot::BezierFitter plot_parabolic_approximation_curve;

    enum class StepKind
    {
      keep_going,               // Just subtract current_learning_rate times derivative.
      done,                     // w was already updated.
      check_energy,             // After adding the new sample, abort if the required energy is too large.
      abort                     // Stop going in the current vdirection.
    };

    StepKind step_kind = StepKind::done;

    while (true)
    {
      // Add new sample to the history.
      Sample const& current = history.add(w, L(w), L.derivative(w), scale);
      current.set_point(plot.create_point(second_layer, point_style, {current.w(), current.Lw()}));
      current.set_label(plot.create_text(second_layer, label_style({.position = cairowindow::draw::centered_below}),
            current.P(), std::to_string(history.total_number_of_samples() - 1)));

      if (step_kind == StepKind::check_energy)
      {
        if (!energy.maybe_update(current.Lw()))
        {
          step_kind = StepKind::abort;
          Dout(dc::notice, "Too much energy used: need to abort this direction.");
        }
      }
      else
      {
        // Update the current kinetic energy. If this is an overshoot, just increase the energy to 0.
        energy.update(current.Lw());

        // The default is that new_sample is set to the extreme of some polynomial.
        step_kind = StepKind::done;
      }

      // Block until a redraw is necessary (for example because the user moved a draggable object,
      // or wants to print the current drawing) then go to the top of loop for a redraw.
      if (!window.handle_input_events())
        break;          // Program must be terminated.

      Dout(dc::notice, "------------------------------------");

      // Erase all previous curves (if they exist).
      plot_approximation_curve.reset();
      plot_derivative_curve.reset();
      plot_quotient_curve.reset();
      plot_parabolic_approximation_curve.reset();

      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      if (step_kind != StepKind::abort)
      {
        //===================================================================================================
        // Create a parabolic approximation from the last two samples (or a line if we only have one sample).
        // If we have only one sample, then new_sample and energy accordingly.

        math::QuadraticPolynomial parabolic_approximation;

        // Do we have only one sample?
        if (history.relevant_samples() == 1)
        {
          double const dLdw = current.dLdw();

          // If we have just one point, then the approximation is a linear function:
          //
          // A(w) = coef[0] + L'(w) w
          parabolic_approximation[1] = dLdw;

          // Did we drop into a (local) minimum as a starting point?!
          if (Scale::almost_zero(current.w(), dLdw))
          {
            // In this case we can't do anything else than just make a step in some random vdirection.
            w -= learning_rate;
          }
          else
          {
            // Just gradient descent: move downhill.
            w -= learning_rate * dLdw;
          }
        }
        else
        {
          // If we have at least two points, the approximation is a parabola:
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
          double inverse_det = 0.5 / (current.w() - prev.w());
          parabolic_approximation[2] = inverse_det * (current.dLdw() - prev.dLdw());
          // See https://math.stackexchange.com/questions/4913175
          parabolic_approximation[1] =
            inverse_det * (2.0 * current.Lw() - 2.0 * prev.Lw() + (prev.dLdw() - current.dLdw()) * (current.w() + prev.w()));
        }
        parabolic_approximation[0] = current.Lw() - parabolic_approximation(current.w());

        // End of determining parabolic approximation.
        //===================================================================================================

#ifdef CWDEBUG
        Dout(dc::notice|continued_cf, "parabolic_approximation = " << parabolic_approximation << " (based on " << current);
        if (history.total_number_of_samples() > 1)
          Dout(dc::continued, " and " << history.prev());
        Dout(dc::finish, ")");
#endif

        // Draw the parabolic approximation.
        plot_parabolic_approximation_curve.solve(
            [&parabolic_approximation](double w) -> Point { return {w, parabolic_approximation(w)}; }, plot.viewport());
        plot.add_bezier_fitter(second_layer, curve_line_style({.line_color = color::red}), plot_parabolic_approximation_curve);

        // Flush all expose events related to the drawing done above.
        window.set_send_expose_events(true);

        //===================================================================================================
        // If we have two or more samples, then update new_sample and energy.

        int relevant_samples = history.relevant_samples();
        if (relevant_samples > 1)
        {
          Sample const& prev = history.prev();
          // β = (L'(w₁) - L'(w₀)) / (w₁ - w₀)    [see README.gradient_descent]
          double beta_inverse = (current.w() - prev.w()) / (current.dLdw() - prev.dLdw());
          if (hdirection == unknown)
          {
            // If beta is negative, then there is no minimum, only a maximum.
            vdirection = beta_inverse < 0.0 ? up : down;
            Dout(dc::notice, "vdirection is set to " << (vdirection == up ? "up" : "down"));
          }
          if (hdirection != unknown && vdirection * beta_inverse < 0.0)
          {
            // There is no extreme in the vdirection that we're going.
            // Just keep going in the same direction as before.
            double step = hdirection * std::abs(*scale);
            if (vdirection == down)
              w -= step;
            else
              w += step;
            Dout(dc::notice, ((vdirection == (step > 0.0 ? up : down)) ? "Incremented" : "Decremented") <<
                " w with scale (" << std::abs(step) << ")");
            // Abort if the result requires more energy than we have.
            step_kind = StepKind::check_energy;
            // Maybe update the scale, if the sample is further away from the vertex and still matches the parabola.
  //FIXME: move to top
//           if (scale->update(new_sample))
//             scale->draw_indicators();
          }
          else
          {
            // Set w to the value where the derivative of this parabolic approximation is zero.
            Dout(dc::notice, "Setting new_sample to the extreme of parabolic approximation:");
            double step = beta_inverse * current.dLdw();
            w -= step;
            Dout(dc::notice, history.current().w() << " --> " << history.total_number_of_samples() << ": " << w);

//FIXME: move to top(?)
            // Initialize or update the scale with the new parabolic approximation.
            if (relevant_samples == 2)
              scale->initialize(prev, current, parabolic_approximation);
            else
            {
              math::QuadraticPolynomial old_parabola = scale->get_parabola();
              if (scale->update(prev, current, parabolic_approximation))
                scale->draw_old_parabola(old_parabola);
            }
            scale->draw_indicators();

            // With the new extreme insight, adjust the learning rate to the new scale.
            current_learning_rate = 0.1 * std::abs(beta_inverse);

            if (relevant_samples > 2)
            {
              double abs_step = std::abs(step);
              Dout(dc::notice, "abs_step = " << abs_step << " (between " <<
                  (history.total_number_of_samples() - 1) << " and " << history.total_number_of_samples() << " (to be added))");

              // Did we reach the (local) extreme?
              if (abs_step < 0.01 * *scale)
              {
                double w2_1 = w;
                double w2_2 = w2_1 * w2_1;
                double w2_3 = w2_2 * w2_1;
                double w2_4 = w2_2 * w2_2;

                // If the new sample is too close to the current one, then ignore current.
                bool skip_sample = std::abs(w2_1 - current.w()) < 0.001 * *scale;
                Sample const& w1 = skip_sample ? prev : current;
                Sample const& w0 = skip_sample ? history.prev_prev() : prev;

                double w1_1 = w1.w();
                double w1_2 = w1_1 * w1_1;
                double w1_3 = w1_2 * w1_1;
                double w1_4 = w1_2 * w1_2;

                double w0_1 = w0.w();
                double w0_2 = w0_1 * w0_1;
                double w0_3 = w0_2 * w0_1;
                double w0_4 = w0_2 * w0_2;

                Dout(dc::notice, "Fitting a fourth degree parabola using the samples at " << w0_1 << ", " << w1_1 << " and " << w2_1);

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
                M <<        1.0,   2.0 *  w2_1,    3.0 * w2_2,    4.0 * w2_3,
                            1.0,   2.0 *  w1_1,    3.0 * w1_2,    4.0 * w1_3,
                            1.0,   2.0 *  w0_1,    3.0 * w0_2,    4.0 * w0_3,
                    w2_1 - w1_1,   w2_2 - w1_2,   w2_3 - w1_3,   w2_4 - w1_4;

                Eigen::Vector4d D;
                //FIXME: do this after adding w to the history?
                D <<                 /*new_sample.dLdw()*/ L.derivative(w),
                                             w1.dLdw(),
                                             w0.dLdw(),
                             /*new_sample.Lw()*/L(w) - w1.Lw();

                Eigen::Vector4d C = M.colPivHouseholderQr().solve(D);

                math::Polynomial approximation(5 COMMA_CWDEBUG_ONLY("w"));
                approximation[1] = C[0];
                approximation[2] = C[1];
                approximation[3] = C[2];
                approximation[4] = C[3];
                approximation[0] = /*new_sample.Lw()*/L(w) - approximation(w);
                Dout(dc::notice, "approximation = " << approximation);

                plot_approximation_curve.solve([&approximation](double w) -> Point { return {w, approximation(w)}; }, plot.viewport());
                plot.add_bezier_fitter(second_layer, curve_line_style({.line_color = color::teal}), plot_approximation_curve);

                auto derivative = approximation.derivative();
                Dout(dc::notice, "derivative = " << derivative);

                plot_derivative_curve.solve([&derivative](double w) -> Point { return {w, derivative(w)}; }, plot.viewport());
                plot.add_bezier_fitter(second_layer, curve_line_style({.line_color = color::magenta}), plot_derivative_curve);

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

                plot_quotient_curve.solve([&quotient](double w) -> Point { return {w, 10.0 * quotient(w)}; }, plot.viewport());
                plot.add_bezier_fitter(second_layer, curve_line_style({.line_color = color::blue}), plot_quotient_curve);

                // Did we reach the (local) extreme?
                if (abs_step < (vdirection == up ? 0.05 : 0.01) * *scale)
                {
                  Dout(dc::notice, (vdirection == up ? "Maximum" : "Minimum") << " reached: " << abs_step << " < 0.01 * " << *scale);

                  // Store the found extreme in the history.
                  //FIXME: can't this be done by top of loop?
                  history.add(w, L(w), L.derivative(w), scale);

                  // Store it as an extreme.
                  std::list<SampleMinimum>::iterator new_extreme;
                  if (hdirection == right)
                    new_extreme = extremes.emplace(extremes.end(), history.current(), scale, energy.energy(), right);
                  else
                    new_extreme = extremes.emplace(extremes.begin(), history.current(), scale, energy.energy(), left);

                  // Keep track of the best minimum so far; or abort if this minimum isn't better then one found before.
                  if (vdirection == down)
                  {
                    if (best_minimum == extremes.end() || best_minimum->Lw() > new_extreme->Lw())
                    {
                      best_minimum = new_extreme;
                      Dout(dc::notice, "best_minimum set to " << best_minimum->w() << " / " << best_minimum->Lw());
                    }
                    if (new_extreme != best_minimum)
                    {
                      // The new minimum isn't better than what we found already. Stop going into this vdirection.
                      step_kind = StepKind::abort;
                    }
                  }

                  // Change vdirection.
                  vdirection = -vdirection;
                  Dout(dc::notice, "vdirection is set to " << (vdirection == up ? "up" : "down") <<
                      "; hdirection = " << (hdirection == left ? "left" : "right"));

                  if (number_of_zeroes > 0)
                  {
                    int best_zero = (number_of_zeroes == 2 &&
                        (hdirection == right || (hdirection == unknown && approximation(zeroes[1]) < approximation(zeroes[0])))) ? 1 : 0;
                    w = zeroes[best_zero];
                    scale->reset(number_of_zeroes < 2 ? *scale : zeroes[1] - zeroes[0]);
                    scale->erase_indicators();
                    history.reset();
                    Dout(dc::notice, "Set w to found extreme: w = " << w);
                    Dout(dc::notice(hdirection == unknown && number_of_zeroes == 2),
                        "Best zero was " << zeroes[best_zero] << " with A(" << zeroes[best_zero] << ") = " <<
                        approximation(zeroes[best_zero]) << " (the other has value " << approximation(zeroes[1 - best_zero]) << ")");

                    // Update the current kinetic energy. If this is an overshoot, just increase the energy to 0.
                  }
                  else if (hdirection == unknown)
                  {
                    if (number_of_zeroes == 0)
                    {
                      // Keep going in the same vdirection.
                      w -= *scale;
  //FIXME: move to top
//           if (scale->update(new_sample))
//             scale->draw_indicators();
                    }

                    // Now that the decision on which hdirection we explore is taken, store that decision.
                    hdirection = w - history.current().w() < 0.0 ? left : right;
                    Dout(dc::notice, "Initializing horizontal direction to " << (hdirection == left ? "left" : "right"));
                  }
                }
              }
            }
          }
        }
      }

      if (step_kind == StepKind::abort)
      {
        // Jump back to the best minimum and continue in the opposite hdirection.

        // Has a minimum been found at all?
        if (best_minimum == extremes.end())
          break;

        // Was this minimum already explored in both directions?
        if (best_minimum->done())
          break;

        // Change hdirection.
        hdirection = -hdirection;
        Dout(dc::notice, "Changed horizontal direction to " << (hdirection == left ? "left" : "right"));

        // Restore the current sample and scale to the values belonging to this minimum.
        w = best_minimum->w();
        scale = best_minimum->scale(hdirection);

        // Forget any acceleration that happened before.
        current_learning_rate = scale->learning_rate();

        // Do a step with a size equal to the scale in this minimum: we expect it to be a parabola up till that point.
        w += hdirection * std::abs(*scale);

        // Restore the energy to what it was when this minimum was stored.
        energy.set(best_minimum->energy(), best_minimum->Lw());
      }

      if (step_kind != StepKind::done)
        Dout(dc::notice, history.current().w() << " --> " << history.total_number_of_samples() << ": " << w);
    }

    if (best_minimum != extremes.end())
      Dout(dc::notice, "Found global minimum " << best_minimum->Lw() << " at w = " << best_minimum->w());
    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
