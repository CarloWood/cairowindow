#include "sys.h"
#include "Polynomial.h"
#include "QuadraticPolynomial.h"
#include "Sample.h"
#include "Scale.h"
#include "Approximation.h"
#include "History.h"
#include "LocalExtreme.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/symbolic/symbolic.h"
#include "utils/square.h"
#include "utils/AIAlert.h"
#include <Eigen/Dense>
#include <sstream>
#include <thread>
#include <vector>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/debug_ostream_operators.h"
#include "utils/print_using.h"
#endif

using utils::has_print_on::operator<<;

#if 1
class Function
{
 public:
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
class PlotSample
{
 private:
  gradient_descent::Sample const* master_{};

  mutable cairowindow::plot::Point P_;
  mutable cairowindow::plot::Text P_label_;

 public:
  // Required for History.
  PlotSample() = default;

  // Constructor used by LocalExtreme.
  PlotSample(gradient_descent::Sample const* master, cairowindow::plot::Point const& point, cairowindow::plot::Text const& label) :
    master_(master), P_(point), P_label_(label) { }

  // Required for History.
  void initialize(gradient_descent::Sample const* master, cairowindow::plot::Point const& point, cairowindow::plot::Text const& label)
  {
    master_ = master;
    P_ = point;
    P_label_ = label;
  }

  cairowindow::Point const& P() const { return P_; }
  cairowindow::Text const& label() const { return P_label_; }

  gradient_descent::Sample const* sample() const { return master_; }

  double w() const { return master_->w(); }
  double Lw() const { return master_->Lw(); }
  double dLdw() const { return master_->dLdw(); }

#ifdef CWDEBUG
  std::string debug_label() const { return P_label_.text(); }

  void print_on(std::ostream& os) const
  {
    os << debug_label() << " (at w = " << master_->w() << ")";
  }
#endif
};

class PlotParabolaScale
{
 public:
  static cairowindow::draw::ConnectorStyle const s_indicator_style;

 private:
  gradient_descent::Scale const& master_;

  // Used to visualize the Scale:
  cairowindow::plot::Plot* plot_ = nullptr;
  boost::intrusive_ptr<cairowindow::Layer> layer_;
  cairowindow::plot::Connector plot_indicator_;
  cairowindow::plot::Text plot_scale_text_;
  cairowindow::plot::Line plot_vertical_line_through_w_;
  cairowindow::plot::Line plot_vertical_line_through_v_;

  // Temporary curves (used while developing this class).
  cairowindow::plot::BezierFitter plot_old_parabola_;

 public:
  PlotParabolaScale(gradient_descent::Scale const& master, cairowindow::plot::Plot& plot,
      boost::intrusive_ptr<cairowindow::Layer> const& layer) :
    master_(master), plot_(&plot), layer_(layer) { }

  void draw_indicators()
  {
    ASSERT(master_.has_sample());

    double x1 = master_.parabola().vertex_x();
    double x2 = master_.edge_sample_w();
    double scale_y = plot_->yrange().min() + 0.5 * plot_->yrange().size();
    plot_indicator_ = cairowindow::plot::Connector{{x1, scale_y}, {x2, scale_y},
        cairowindow::Connector::open_arrow, cairowindow::Connector::open_arrow};
    plot_->add_connector(layer_, s_indicator_style, plot_indicator_);
    plot_scale_text_ = plot_->create_text(layer_, {{.position = cairowindow::draw::centered_above, .offset = 2.0}},
        cairowindow::Point{(x1 + x2) / 2, scale_y}, "scale");

    plot_vertical_line_through_v_ = cairowindow::plot::Line{{master_.parabola().vertex_x(), scale_y}, cairowindow::Direction::up};
    plot_->add_line(layer_, s_indicator_style({.dashes = {3.0, 3.0}}), plot_vertical_line_through_v_);

    plot_vertical_line_through_w_ = cairowindow::plot::Line{{master_.edge_sample_w(), scale_y}, cairowindow::Direction::up};
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
    return master_.parabola();
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

class PlotHistory : public gradient_descent::History
{
 private:
  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> const& layer_;
  cairowindow::draw::PointStyle const& point_style_;
  cairowindow::draw::TextStyle const& label_style_;
  utils::Array<PlotSample, size, gradient_descent::HistoryIndex> plot_samples_;

 public:
  PlotHistory(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer,
      cairowindow::draw::PointStyle const& point_style, cairowindow::draw::TextStyle const& label_style) :
    plot_(plot), layer_(layer), point_style_(point_style), label_style_(label_style) { }

  PlotSample const& add(double w, double Lw, double dLdw, gradient_descent::Scale const& scale, bool& current_is_replacement)
  {
    gradient_descent::HistoryIndex index = gradient_descent::History::add(w, Lw, dLdw, scale, current_is_replacement);
    plot_samples_[index].initialize(&current(),
      plot_.create_point(layer_, point_style_, {w, Lw}),
      plot_.create_text(layer_, label_style_({.position = cairowindow::draw::centered_below}),
            cairowindow::Point{w, Lw}, std::to_string(total_number_of_samples() - 1)));
    return plot_samples_[index];
  }

  cairowindow::plot::Plot& plot() const { return plot_; }
  boost::intrusive_ptr<cairowindow::Layer> const& layer() const { return layer_; }
  cairowindow::draw::PointStyle const& point_style() const { return point_style_; }
  cairowindow::draw::TextStyle const& label_style() const { return label_style_; }
};

//static
cairowindow::draw::ConnectorStyle const PlotParabolaScale::s_indicator_style{{.line_width = 1}};

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

struct PlotLocalExtreme : gradient_descent::LocalExtreme
{
  PlotParabolaScale plot_approximation_parabola_scale_;         // Sibling of gradient_descent::LocalExtreme::approximation_
  PlotSample plot_vertex_sample_;                               // Sibling of gradient_descent::LocalExtreme::vertex_sample_

  PlotLocalExtreme(PlotHistory const& history, gradient_descent::Approximation const& approximation, double energy) :
    gradient_descent::LocalExtreme(history.current(), approximation, energy),
    plot_approximation_parabola_scale_(approximation_.parabola_scale(), history.plot(), history.layer()),
    plot_vertex_sample_(&vertex_sample_,
        history.plot().create_point(history.layer(), history.point_style()({.color_index = 15}), {vertex_sample_.w(), vertex_sample_.Lw()}),
        history.plot().create_text(history.layer(), history.label_style()({.position = cairowindow::draw::centered_below}),
          cairowindow::Point{vertex_sample_.w(), vertex_sample_.Lw()}, std::to_string(history.total_number_of_samples() - 1))) { }
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
    using namespace gradient_descent;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Gradient descent of " + L.as_string(), 600, 450);

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
    PlotHistory history(plot, second_layer, point_style, label_style);

    // Initial value.
    Weight w(w_0);
    Dout(dc::notice, "Initial value of new_sample; w = " << w);

    Approximation current_approximation;
    PlotParabolaScale current_plot_approximation_parabola_scale(current_approximation.parabola_scale(), plot, second_layer);
    Approximation* approximation_ptr = &current_approximation;
    PlotParabolaScale* plot_approximation_parabola_scale_ptr = &current_plot_approximation_parabola_scale;

    KineticEnergy energy(plot, second_layer, L_max);
    Dout(dc::notice, "Initial height is " << energy);
    constexpr int unknown = 0;
    constexpr int down = 1;
    constexpr int up = -1;
    int vdirection = down;
    HorizontalDirection hdirection = unknown_horizontal_direction;
    // Pointers to elements may not be invalidated.
    using extremes_type = std::list<PlotLocalExtreme>;
    extremes_type extremes;
    extremes_type::iterator best_minimum = extremes.end();
    extremes_type::iterator last_extreme = extremes.end();

    plot::BezierFitter plot_approximation_curve;
    plot::BezierFitter plot_derivative_curve;
    plot::BezierFitter plot_quotient_curve;
    plot::BezierFitter plot_fourth_degree_approximation_curve;

    enum class StepKind
    {
      done,                     // w was already updated.
      check_energy,             // After adding the new sample, abort if the required energy is too large.
      abort                     // Stop going in the current vdirection.
    };

    StepKind step_kind = StepKind::done;

    while (true)
    {
      // Add new sample to the history.
      bool current_is_replacement;
      PlotSample const& current = history.add(w, L(w), L.derivative(w), approximation_ptr->parabola_scale(), current_is_replacement);
      // If current_is_replacement is true then the previous `current` had a value so close to w that
      // it was replaced instead of adding a new sample to the history.

      if (step_kind == StepKind::check_energy)
      {
        if (!energy.maybe_update(current.Lw()))
        {
          step_kind = StepKind::abort;
          Dout(dc::notice, "Too much energy used: need to abort this direction.");
          ASSERT(hdirection != unknown_horizontal_direction);
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
      plot_fourth_degree_approximation_curve.reset();

      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

      if (step_kind != StepKind::abort)
      {
        //===================================================================================================
        // Create a parabolic approximation from the last two samples (or a line if we only have one sample).
        // If we have only one sample, then new_sample and energy accordingly.

        Approximation& approximation(*approximation_ptr);
        PlotParabolaScale& plot_parabola_scale(*plot_approximation_parabola_scale_ptr);

        ScaleUpdate result = approximation.add(current.sample(), current_is_replacement);
        math::QuadraticPolynomial old_parabola = plot_parabola_scale.get_parabola();

        switch (result)
        {
          case ScaleUpdate::first_sample:
            break;
          case ScaleUpdate::initialized:
            plot_parabola_scale.draw_indicators();
            break;
          case ScaleUpdate::towards_vertex:
            plot_parabola_scale.draw_indicators();
            plot_parabola_scale.draw_old_parabola(old_parabola);
            break;
          case ScaleUpdate::away_from_vertex:
            plot_parabola_scale.draw_indicators();
            break;
        }

        // Do we have only one sample?
        if (history.relevant_samples() == 1)
        {
          double const dLdw = current.dLdw();

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

        // End of determining parabolic approximation.
        //===================================================================================================

#ifdef CWDEBUG
        Dout(dc::notice, "approximation = " << approximation <<
            " (" << utils::print_using(approximation, &Approximation::print_based_on) << ")");
#endif

        // Draw the parabolic approximation.
        plot_approximation_curve.solve(
            [&approximation](double w) -> Point { return {w, approximation.parabola()(w)}; }, plot.viewport());
        plot.add_bezier_fitter(second_layer, curve_line_style({.line_color = color::red}), plot_approximation_curve);

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
          if (hdirection == unknown_horizontal_direction)
          {
            // If beta is negative, then there is no minimum, only a maximum.
            vdirection = beta_inverse < 0.0 ? up : down;
            Dout(dc::notice, "vdirection is set to " << (vdirection == up ? "up" : "down"));
          }
          if (hdirection != unknown_horizontal_direction && vdirection * beta_inverse < 0.0)
          {
            // There is no extreme in the vdirection that we're going.
            // Just keep going in the same direction as before.
            double step = static_cast<int>(hdirection) * std::abs(approximation.parabola_scale());
            if (vdirection == down)
              w -= step;
            else
              w += step;
            Dout(dc::notice, ((vdirection == (step > 0.0 ? up : down)) ? "Incremented" : "Decremented") <<
                " w with scale (" << std::abs(step) << ")");
            // Abort if the result requires more energy than we have.
            step_kind = StepKind::check_energy;
          }
          else
          {
            // Set w to the value where the derivative of this parabolic approximation is zero.
            Dout(dc::notice, "Setting new_sample to the extreme of parabolic approximation:");
            double step = w;
            w = approximation.parabola().vertex_x();
            step -= w;
#if 0
            // Use the actual derivative of L(w), instead of the derivative of the parabolic approximation
            // which might be a little different when not both points lay on a piece of the curve that
            // actually is a parabola.
            //
            // As a result, this new w might not be exactly equal to the vertex of the approximated parabola.
            double step = beta_inverse * current.dLdw();
            w -= step;
#endif
            Dout(dc::notice, history.current().w() << " --> " << history.total_number_of_samples() << ": " << w);

            if (relevant_samples > 2)
            {
              double abs_step = std::abs(step);
              Dout(dc::notice, "abs_step = " << abs_step << " (between " <<
                  (history.total_number_of_samples() - 1) << " and " << history.total_number_of_samples() << " (to be added))");

              // Did we reach the (local) extreme?
              if (abs_step < 0.01 * approximation.parabola_scale())
              {
                double w2_1 = w;
                double w2_2 = w2_1 * w2_1;
                double w2_3 = w2_2 * w2_1;
                double w2_4 = w2_2 * w2_2;

                // If the new sample is too close to the current one, then ignore current.
                bool skip_sample = std::abs(w2_1 - current.w()) < 0.001 * approximation.parabola_scale();
                Sample const& w1 = skip_sample ? prev : *current.sample();
                Sample const& w0 = skip_sample ? history.prev_prev() : prev;

                double w1_1 = w1.w();
                double w1_2 = w1_1 * w1_1;
                double w1_3 = w1_2 * w1_1;
                double w1_4 = w1_2 * w1_2;

                double w0_1 = w0.w();
                double w0_2 = w0_1 * w0_1;
                double w0_3 = w0_2 * w0_1;
                double w0_4 = w0_2 * w0_2;

                Dout(dc::notice, "Fitting a fourth degree polynomial using the samples at " << w0_1 << ", " << w1_1 << " and " << w2_1);

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

                math::Polynomial fourth_degree_approximation(5 COMMA_CWDEBUG_ONLY("w"));
                fourth_degree_approximation[1] = C[0];
                fourth_degree_approximation[2] = C[1];
                fourth_degree_approximation[3] = C[2];
                fourth_degree_approximation[4] = C[3];
                fourth_degree_approximation[0] = /*new_sample.Lw()*/L(w) - fourth_degree_approximation(w);
                Dout(dc::notice, "approximation = " << fourth_degree_approximation);

                plot_fourth_degree_approximation_curve.solve(
                    [&fourth_degree_approximation](double w) -> Point { return {w, fourth_degree_approximation(w)}; }, plot.viewport());
                plot.add_bezier_fitter(second_layer, curve_line_style({.line_color = color::teal}), plot_fourth_degree_approximation_curve);

                auto derivative = fourth_degree_approximation.derivative();
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

                // Isn't this always true?!
                ASSERT(abs_step < (vdirection == up ? 0.05 : 0.01) * approximation.parabola_scale());

                // Did we reach the (local) extreme?
                if (abs_step < (vdirection == up ? 0.05 : 0.01) * approximation.parabola_scale())
                {
                  Dout(dc::notice, (vdirection == up ? "Maximum" : "Minimum") << " reached: " << abs_step <<
                      " < 0.01 * " << approximation.parabola_scale());

                  // Store the found extreme in the history.
                  //FIXME: can't this be done by top of loop?
                  history.add(w, L(w), L.derivative(w), approximation.parabola_scale(), current_is_replacement);

                  // Following adding the first extreme, we need to decide on a horizontal direction!
                  ASSERT(hdirection != unknown_horizontal_direction || extremes.empty());

                  // Store it as an extreme.
                  extremes_type::iterator new_extreme;
                  if (hdirection == right)
                    new_extreme = extremes.emplace(extremes.end(), history, approximation, energy.energy());
                  else
                    new_extreme = extremes.emplace(extremes.begin(), history, approximation, energy.energy());

                  // Switch approximation_ptr to the parabolic approximation stored in this extreme:
                  // we need to keep updating it when new samples are added that match the same parabolic.
                  plot_approximation_parabola_scale_ptr->erase_indicators();
                  approximation_ptr = &new_extreme->approximation(hdirection);
                  plot_approximation_parabola_scale_ptr = &new_extreme->plot_approximation_parabola_scale_;

                  // Keep track of the best minimum so far; or abort if this minimum isn't better then one found before.
                  if (vdirection == down)
                  {
                    if (best_minimum == extremes.end() || best_minimum->vertex_sample().Lw() > new_extreme->vertex_sample().Lw())
                    {
                      best_minimum = new_extreme;
                      Dout(dc::notice, "best_minimum set to " << best_minimum->vertex_sample() <<
                          " and parabolic approximation: " << approximation);
                    }
                    if (new_extreme != best_minimum)
                    {
                      // The new minimum isn't better than what we found already. Stop going into this direction.
                      Dout(dc::notice, "The new minimum isn't better than what we found already. Stop going into the direction " << hdirection);
                      step_kind = StepKind::abort;
                    }
                  }

                  ASSERT(extremes.size() != 1 ||
                      (vdirection == down && best_minimum != extremes.end() && step_kind != StepKind::abort));

                  // Change vdirection.
                  vdirection = -vdirection;
                  Dout(dc::notice, "vdirection is set to " << (vdirection == up ? "up" : "down") << " (hdirection = " << hdirection << ")");

                  if (number_of_zeroes > 0)
                  {
                    int best_zero = (number_of_zeroes == 2 &&
                        (hdirection == right ||
                         (hdirection == unknown_horizontal_direction &&
                          fourth_degree_approximation(zeroes[1]) < fourth_degree_approximation(zeroes[0])))) ? 1 : 0;
                    w = zeroes[best_zero];
                    // Use the Approximation object on the stack again (as opposed to one from a LocalExtreme).
                    plot_approximation_parabola_scale_ptr->erase_indicators();
                    approximation_ptr = &current_approximation;
                    plot_approximation_parabola_scale_ptr = &current_plot_approximation_parabola_scale;
                    // Reset the parabolic approximation.
                    approximation_ptr->reset();
                    history.reset();
                    Dout(dc::notice, "Set w to found extreme: w = " << w);
                    Dout(dc::notice(hdirection == unknown_horizontal_direction && number_of_zeroes == 2),
                        "Best zero was " << zeroes[best_zero] << " with A(" << zeroes[best_zero] << ") = " <<
                        fourth_degree_approximation(zeroes[best_zero]) <<
                        " (the other has value " << fourth_degree_approximation(zeroes[1 - best_zero]) << ")");

                    // Update the current kinetic energy. If this is an overshoot, just increase the energy to 0.
                  }

                  if (hdirection == unknown_horizontal_direction)
                  {
                    if (number_of_zeroes == 0)
                    {
                      // The "scale" of the (current) parabolic approximation is set to 'edge sample' minus x-coordinate of the vertex (v_x),
                      // where 'edge sample' is the sample that is "part of" the parabolic approximation that is the furthest away
                      // from the vertex. "Part of" here means that it deviates vertically less than 10% of the vertical distance to
                      // the vertex. In that sense we can consider that we "came from" the direction of that edge sample.
                      // For example, if the edge sample had x-coordinate w = w_E and we came from the right (w_E > v_x) then
                      // scale = w_E - v_x > 0. Subtracting a positive scale thus means we continue in the direction towards to the left.
                      // If w_E is on the left of the vertex then scale is negative and subtracting is causes us to continue to the right.
                      //
                      // Keep going in the same vdirection.
                      w -= approximation.parabola_scale();

                      // Note that w was already set to the v_x before, but the test below still works
                      // because |v_x - history.current().w()| was determined to be less than 1% of approximation.parabola_scale().
                    }

                    // Now that the decision on which hdirection we explore is taken, store that decision.
                    hdirection = w - history.current().w() < 0.0 ? left : right;
                    Dout(dc::notice, "Initializing horizontal direction to " << hdirection);
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

        ASSERT(hdirection != unknown_horizontal_direction);

        // Change hdirection.
        hdirection = opposite(hdirection);
        Dout(dc::notice, "Changed horizontal direction to " << hdirection);

        // Restore the current sample and scale to the values belonging to this minimum.
        w = best_minimum->vertex_sample().w();
        Dout(dc::notice, "Restored w to best minimum at " << w);

        // Do a step with a size equal to the scale in this minimum: we expect it not to drastically change before that point.
        w += static_cast<int>(hdirection) * std::abs(best_minimum->approximation(hdirection).parabola_scale());

        // Restore the energy to what it was when this minimum was stored.
        energy.set(best_minimum->energy(), best_minimum->vertex_sample().Lw());

        // Use the Approximation object on the stack again (as opposed to one from a LocalExtreme).
        plot_approximation_parabola_scale_ptr->erase_indicators();
        approximation_ptr = &current_approximation;
        plot_approximation_parabola_scale_ptr = &current_plot_approximation_parabola_scale;
        // Reset the parabolic approximation.
        approximation_ptr->reset();
      }

      if (step_kind != StepKind::done)
        Dout(dc::notice, history.current().w() << " --> " << history.total_number_of_samples() << ": " << w);
    }

    if (best_minimum != extremes.end())
      Dout(dc::notice, "Found global minimum " << best_minimum->vertex_sample());
    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
