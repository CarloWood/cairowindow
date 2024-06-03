#include "sys.h"
#include "Polynomial.h"
#include "QuadraticPolynomial.h"
#include "Sample.h"
#include "Scale.h"
#include "Approximation.h"
#include "History.h"
#include "LocalExtreme.h"
#include "KineticEnergy.h"
#include "HorizontalDirection.h"
#include "VerticalDirection.h"
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

#define USE_SLIDERS 0

using utils::has_print_on::operator<<;

class FunctionBase
{
 protected:
  symbolic::Symbol const& w_ = symbolic::Symbol::realize("w");
  symbolic::Symbol const& a_ = symbolic::Symbol::realize("a");
  symbolic::Symbol const& b_ = symbolic::Symbol::realize("b");
  symbolic::Symbol const& c_ = symbolic::Symbol::realize("c");
  symbolic::Symbol const& d_ = symbolic::Symbol::realize("d");
  symbolic::Symbol const& e_ = symbolic::Symbol::realize("e");
};

#if 1
class Function : FunctionBase
{
 public:
  static constexpr double w_0 = -51.0;
  static constexpr double w_min = -60.0; //-80.0;
  static constexpr double w_max = -50.0; //20.0;

 private:
  using Expression = symbolic::Expression;
  using Constant = symbolic::Constant;
  using Symbol = symbolic::Symbol;
  using Func = symbolic::Function;

  static Constant const tp;
  Symbol const& amplitude_ = symbolic::Symbol::realize("amplitude");
  Symbol const& level_ = symbolic::Symbol::realize("level");
  Symbol const& phase_ = symbolic::Symbol::realize("phase");

 private:
  Func const& sigmoid_ = Func::realize("sigmoid", exp(3 * (w_ + tp)) / (1 + exp(3 * (w_ + tp))));

  Expression const& function_ = (1.0 - sigmoid_) * (a_ + b_ * w_ + c_ * (w_^2)) +
    (sigmoid_ * (amplitude_ * exp((tp - w_) / 10) * sin(d_ * w_ + phase_) + level_));

  Expression const& derivative_ = function_.derivative(w_);

 public:
  Function()
  {
    a_ = 14.6;
    b_ = 3.15;
    c_ = 0.451;
    d_ = 0.5;
    amplitude_ = 0.012027;
    level_ = 1878.38;
    phase_ = 1.91892;
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
    return function_.evaluate();
  }

  double derivative(double w) const
  {
    w_ = w;
    sigmoid_.reset_evaluation();
    return derivative_.evaluate();
  }

#if USE_SLIDERS
 public:
  void set_sliders(double amplitude, double level, double phase)
  {
    amplitude_ = amplitude;
    level_ = level;
    phase_ = phase;
    sigmoid_.reset_evaluation();
  }
#endif
};

//static
symbolic::Constant const Function::tp = symbolic::Constant::realize(55);

#elif 0
class Function : FunctionBase
{
 public:
  static constexpr double w_0 = 12.0;
  static constexpr double w_min = -20.0;
  static constexpr double w_max = 80.0;

 private:
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
    a2_ = 14.6;
    b2_ = 3.15;
    c2_ = 0.451;
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

  gradient_descent::Scale const& scale() const
  {
    return master_;
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

  void add(double w, double Lw, double dLdw, gradient_descent::Scale const& scale, bool& current_is_replacement)
  {
    gradient_descent::HistoryIndex index = gradient_descent::History::add(w, Lw, dLdw, scale, current_is_replacement);
    plot_samples_[index].initialize(&current(),
      plot_.create_point(layer_, point_style_, {w, Lw}),
      plot_.create_text(layer_, label_style_({.position = cairowindow::draw::centered_below}),
            cairowindow::Point{w, Lw}, std::to_string(total_number_of_samples() - 1)));
  }

  cairowindow::plot::Plot& plot() const { return plot_; }
  boost::intrusive_ptr<cairowindow::Layer> const& layer() const { return layer_; }
  cairowindow::draw::PointStyle const& point_style() const { return point_style_; }
  cairowindow::draw::TextStyle const& label_style() const { return label_style_; }
};

//static
cairowindow::draw::ConnectorStyle const PlotParabolaScale::s_indicator_style{{.line_width = 1}};

class PlotKineticEnergy : public gradient_descent::KineticEnergy
{
  static cairowindow::draw::ConnectorStyle const s_indicator_style;

 private:
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

 public:
  PlotKineticEnergy(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer, double Lw) :
    gradient_descent::KineticEnergy(Lw), plot_(plot), layer_(layer) { }

  void set(double max_Lw, double Lw)
  {
    gradient_descent::KineticEnergy::set(max_Lw, Lw);
    draw_indicators();
  }

  bool maybe_update(double new_Lw)
  {
    if (!gradient_descent::KineticEnergy::maybe_update(new_Lw))
      return false;
    draw_indicators();
    return true;
  }

  void update(double new_Lw)
  {
    gradient_descent::KineticEnergy::update(new_Lw);
    draw_indicators();
  }
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
cairowindow::draw::ConnectorStyle const PlotKineticEnergy::s_indicator_style{{.line_width = 1}};

class AcceleratedGradientDescent
{
  using Approximation = gradient_descent::Approximation;
  using HorizontalDirection = gradient_descent::HorizontalDirection;
  using VerticalDirection = gradient_descent::VerticalDirection;
  using Weight = gradient_descent::Weight;
  using Sample = gradient_descent::Sample;
  using Scale = gradient_descent::Scale;
  using HistoryIndex = gradient_descent::HistoryIndex;

  static constexpr cairowindow::draw::LineStyle curve_line_style_{{.line_width = 1.0}};
  static cairowindow::draw::ConnectorStyle s_difference_expected_style;

 private:
  double learning_rate_;        // In unit_of(w)^2 / unit_of(L).
  double small_step_{};         // This will replace learning_rate_ as soon as we have an idea of the scale of changes.

  // Remember the (most recent) history of samples.
  PlotHistory history_;
  HistoryIndex clamped_history_index_;
  Approximation current_approximation_;
  PlotParabolaScale current_plot_approximation_parabola_scale_;
  Approximation* approximation_ptr_;
  PlotParabolaScale* plot_approximation_parabola_scale_ptr_;
  PlotKineticEnergy energy_;
  double expected_Lw_;                  // Whenever w is changed, this is set to what Lw value the approximation is expecting there.

  // hdirection_ is set when we find a local minimum and decide to explore left or right of that.
  // The result is that we'll rather go away from the vertex of the current matching parabolic approximation
  // then towards it, if that doesn't match the current hdirection_.
  HorizontalDirection hdirection_;
  VerticalDirection vdirection_;

  // Using a std::list: pointers to elements may never be invalidated.
  using extremes_type = std::list<PlotLocalExtreme>;
  extremes_type extremes_;
  extremes_type::iterator best_minimum_;
  extremes_type::iterator last_extreme_;

  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> const& layer_;
  cairowindow::plot::BezierFitter plot_approximation_curve_;
  cairowindow::plot::BezierFitter plot_derivative_curve_;
  cairowindow::plot::BezierFitter plot_quotient_curve_;
  cairowindow::plot::BezierFitter plot_fourth_degree_approximation_curve_;
  cairowindow::plot::Connector plot_difference_expected_;

  enum class IterationState
  {
    done,                     // w was already updated.
    check_energy,             // After adding the new sample, abort if the required energy is too large.
    vertex_jump,              // Only add the new sample to the history if we didn't overshoot another extreme dramatically.
    clamped,                  // Attemped to jump to a vertex, but against hdirection, passed a previous sample. Use that last sample instead.
    gamma_based,              // Internal state used to signify that w was already updated and a call to handle_parabolic_approximation
                              // is not longer desired.
    local_extreme,            // After adding the new sample, handle the fact that we found an extreme.
    abort_hdirection          // Stop going in the current hdirection_.
  };

  IterationState state_{IterationState::done};

 public:
  AcceleratedGradientDescent(double learning_rate, double L_max,
      cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer,
      cairowindow::draw::PointStyle const& point_style, cairowindow::draw::TextStyle const& label_style) :
    learning_rate_(learning_rate),
    history_(plot, layer, point_style, label_style),
    current_plot_approximation_parabola_scale_(current_approximation_.parabola_scale(), plot, layer),
    approximation_ptr_(&current_approximation_),
    plot_approximation_parabola_scale_ptr_(&current_plot_approximation_parabola_scale_),
    energy_(plot, layer, L_max),
    best_minimum_(extremes_.end()),
    last_extreme_(extremes_.end()),
    hdirection_(HorizontalDirection::undecided),
    vdirection_(VerticalDirection::down),
    plot_(plot),
    layer_(layer)
  {
  }

  bool operator()(Weight& w, double Lw, double dLdw);
  bool update_energy();
  void reset_history();
  bool handle_local_extreme(Weight& w);
  void update_approximation(bool current_is_replacement);
  void handle_single_sample(Weight& w);
  void handle_parabolic_approximation(Weight& w);
  bool handle_abort_hdirection(Weight& w);

  bool success() const
  {
    return best_minimum_ != extremes_.end();
  }

  Sample const& minimum() const
  {
    // Only call this function when success() returns true.
    ASSERT(success());
    return best_minimum_->vertex_sample();
  }
};

//static
cairowindow::draw::ConnectorStyle AcceleratedGradientDescent::s_difference_expected_style{{.line_color = cairowindow::color::blue,
  .line_width = 1.0}};

bool AcceleratedGradientDescent::operator()(Weight& w, double Lw, double dLdw)
{
  DoutEntering(dc::notice, "AcceleratedGradientDescent::operator()(" << w << ", " << Lw << ", " << dLdw << ")");

  using namespace gradient_descent;

  double gamma_based_w;
  if (state_ == IterationState::vertex_jump)
  {
    Approximation& approximation(*approximation_ptr_);
    // How else could we have made a parabolic approximation?
    ASSERT(approximation.number_of_relevant_samples() > 0);
    // Then state should have been IterationState::local_extreme.
    ASSERT(!approximation.is_extreme());

    // See https://math.stackexchange.com/questions/4923841/ plus answer.

    double h = Lw - expected_Lw_;
    double c = approximation.parabola()[2];
    double w0 = approximation.current().w();
    double w1 = w;

    double gamma = 6.0 * h / (c * utils::square(w0 - w1));

    Dout(dc::notice, "h = " << h << "; gamma = " << gamma);

    // Did we substantially miss our target?
    if (gamma > -1.0 && std::abs(gamma) > 0.1)
    {
      // Approximate a better w based on the minimum of a third degree polynomial fit, using the new sample.
      gamma_based_w = w0 + 2.0 * (std::sqrt(1.0 + gamma) - 1.0) / gamma * (w1 - w0);
      //expected_Lw_ = ; Not set... is currently ignored anyway.
      if (gamma > 1.207)
      {
        // Too far off, don't even add w to the history.
        Dout(dc::notice, w << " --> " << history_.total_number_of_samples() << ": " <<
            gamma_based_w << " [gamma based] [expected_Lw: " << expected_Lw_ << "]");
        w = gamma_based_w;
        state_ = IterationState::done;
        Dout(dc::notice, "Returning: " << w);
        return true;
      }
      Dout(dc::notice, w << " --> " << (history_.total_number_of_samples() + 1) << ": " <<
          gamma_based_w << " [gamma based] [expected_Lw: " << expected_Lw_ << "]");
      state_ = IterationState::gamma_based;
    }
  }

  // If the new sample (w) is too close to the previous sample (Scale::negligible returns true)
  // then the new sample replaces the previous sample, current_is_replacement is set to true
  // and no new sample is added to the history.
  //
  // Otherwise, the new sample is added to the history and current_is_replacement is set to false.
  bool current_is_replacement;
  history_.add(w, Lw, dLdw, approximation_ptr_->parabola_scale(), current_is_replacement);

  // This function should never return a value whose difference with the previous sample is negligible
  // if there is only single relevant sample in the history.
  ASSERT(!current_is_replacement || history_.relevant_samples() > 1);

  // Erase all previous curves (if they exist).
  plot_approximation_curve_.reset();
  plot_derivative_curve_.reset();
  plot_quotient_curve_.reset();
  plot_fourth_degree_approximation_curve_.reset();

  // Plot the vertical difference from what we expected to what we got.
  plot_difference_expected_ = cairowindow::plot::Connector{{w, expected_Lw_}, {w, Lw},
      cairowindow::Connector::no_arrow, cairowindow::Connector::open_arrow};
  plot_.add_connector(layer_, s_difference_expected_style, plot_difference_expected_);

  // Update kinetic energy. Returns false if too much energy was used.
  if (!update_energy())
    return handle_abort_hdirection(w);

  // Handle the case where this sample is a local extreme.
  if (state_ == IterationState::local_extreme)
  {
    // This state is set when locally the curve looks like a parabola that we're trying to find the
    // extreme of, and the last sample is that extreme (that is, the derivative is now close to zero).
    // Returns false if this exreme is a minimum but isn't better than the previously found best minimum.
    if (!handle_local_extreme(w))
      return handle_abort_hdirection(w);
    Dout(dc::notice, "Returning: " << w);
    return true;        // w was successfully updated by handle_local_extreme.
  }

  // Create/update a parabolic approximation from this and the previous sample (or a line if this is the first sample).
  update_approximation(current_is_replacement);
  Approximation& approximation(*approximation_ptr_);

  if (approximation.number_of_relevant_samples() == 1)
  {
    handle_single_sample(w);
    Dout(dc::notice, "Returning: " << w);
    return true;
  }

  if (hdirection_ == HorizontalDirection::undecided)
  {
    // Do a step towards the extreme of the parabolic approximation.
    vdirection_ = approximation.has_maximum() ? VerticalDirection::up : VerticalDirection::down;
    Dout(dc::notice, "vdirection_ is set to " << vdirection_ << " because hdirection_ is still undecided.");
  }

  if (state_ == IterationState::gamma_based)
  {
    w = gamma_based_w;
    state_ = IterationState::done;
  }
  else
  {
    handle_parabolic_approximation(w);
    if (state_ == IterationState::clamped)
    {
      Dout(dc::notice, "Clamped by sample " << history_[clamped_history_index_]);
      // The vertex of the current approximation doesn't make sense, aka the approximation doesn't make sense.
      // We need to find a reasonable jump point based on the current sample and the one given by clamped_history_index_.
      //FIXME: make this 0.5
      w = 0.52 * (w + history_[clamped_history_index_].w());
      state_ = IterationState::done;
    }
  }

  Dout(dc::notice, "Returning: " << w);
  return true;
}

bool AcceleratedGradientDescent::update_energy()
{
  // Get the height reached.
  double Lw = history_.current().Lw();

  if (state_ == IterationState::check_energy)
  {
    // Update the current kinetic energy. If this is an overshoot, abort this horizontal direction.
    //
    // This state is set when locally the curve looks like a parabola with its minimum
    // on the side that we just came from: in this case we move away from the vertex,
    // losing energy. Since we're going uphill, we need to check if we didn't go too
    // high for the kinetic energy that we have, in which case the exploration of this
    // direction is aborted.

    // If the horizontal direction is still unknown, then we should always go towards the extreme
    // of the local parabolic approximation: in fact we would have jumped there and would not
    // care about the amount of available energy.
    ASSERT(hdirection_ != HorizontalDirection::undecided);

    // Update the energy and check if we had enough energy to reach this height.
    if (!energy_.maybe_update(Lw))
    {
      Dout(dc::notice, "Too much energy used: need to abort this direction (" << hdirection_ << ").");
      return false;
    }
  }
  else
  {
    // Update the current kinetic energy. If this is an overshoot, just increase the energy to 0.
    energy_.update(Lw);
  }
  return true;
}

void AcceleratedGradientDescent::reset_history()
{
  DoutEntering(dc::notice, "AcceleratedGradientDescent::reset_history()");

  // Use the Approximation object on the stack again (as opposed to one from a LocalExtreme).
  plot_approximation_parabola_scale_ptr_->erase_indicators();
  approximation_ptr_ = &current_approximation_;
  plot_approximation_parabola_scale_ptr_ = &current_plot_approximation_parabola_scale_;
  // Reset the parabolic approximation.
  approximation_ptr_->reset();
  history_.reset();
}

bool AcceleratedGradientDescent::handle_local_extreme(Weight& w)
{
  // Following adding the first extreme, we should have decided on a horizontal direction.
  ASSERT(hdirection_ != HorizontalDirection::undecided || extremes_.empty());

  // Store it as an extreme.
  extremes_type::iterator new_extreme =
    extremes_.emplace(hdirection_ == HorizontalDirection::right ? extremes_.end() : extremes_.begin(),
        history_, *approximation_ptr_, energy_.energy());

  // Switch approximation_ptr to the parabolic approximation stored in this extreme:
  // we need to keep updating it when new samples are added that match the same parabolic.
  plot_approximation_parabola_scale_ptr_->erase_indicators();
  approximation_ptr_ = &new_extreme->approximation();
  plot_approximation_parabola_scale_ptr_ = &new_extreme->plot_approximation_parabola_scale_;

  // If we came from (say) the left, and already found a minimum there, then mark left as explored.
  if (best_minimum_ != extremes_.end())
    new_extreme->explored(opposite(hdirection_));

  // Update small_step_ (parabola_scale() returns an absolute value).
  small_step_ = approximation_ptr_->parabola_scale();
  Dout(dc::notice, "small_step_ set to " << small_step_);

  // Keep track of the best minimum so far; or abort if this minimum isn't better then one found before.
  if (vdirection_ == VerticalDirection::down)
  {
    // With vdirection_ down, we were looking for a minimum.
    ASSERT(new_extreme->is_minimum());
    if (best_minimum_ == extremes_.end() || best_minimum_->vertex_sample().Lw() > new_extreme->vertex_sample().Lw())
    {
      best_minimum_ = new_extreme;
      Dout(dc::notice, "best_minimum_ set to " << best_minimum_->vertex_sample() <<
          " and parabolic approximation: " << best_minimum_->approximation());
    }
    if (new_extreme != best_minimum_)
    {
      // The new minimum isn't better than what we found already. Stop going into this direction.
      Dout(dc::notice, "The new minimum (at " << new_extreme->vertex_sample() << ") isn't better than what we found already. "
          "Stop going into the direction " << hdirection_ << ".");
      state_ = IterationState::abort_hdirection;
      return false;
    }
  }

  // After finding a maximum we want to find a minimum and visa versa. Change vdirection_.
  vdirection_ = opposite(vdirection_);
  Dout(dc::notice, "vdirection_ is toggled to " << vdirection_ << ".");

  Sample const& w2 = history_.current();
  double w2_1 = w2.w();

  // Find all samples that are within the scale range of the current approximation.
  double const scale = approximation_ptr_->parabola_scale();
  std::array<int, 5> usable_samples;
  int number_of_usable_samples = 0;
  for (int i = 1; i < history_.relevant_samples() && number_of_usable_samples < usable_samples.size(); ++i)
  {
    Sample const& sample = history_.prev(i);
    double dist = std::abs(sample.w() - w2_1);
    // If the current sample is too close or too far away from the vertex, then skip this sample.
    if (dist < 0.001 * scale || dist > 1.1 * scale)
      continue;

    usable_samples[number_of_usable_samples++] = i;
  }
  Dout(dc::notice, "Number of samples within scale range: " << number_of_usable_samples);
  // If not enough samples we need to get another one! (to be implemented)
  ASSERT(number_of_usable_samples >= 2);
  // Brute force find the two samples that, together with the current sample, have the largest spread.
  int i0 = 0;
  int i1 = 1;
  double best_spread = 0.0;
  if (number_of_usable_samples > 2)
  {
    for (int t0 = 0; t0 < number_of_usable_samples - 1; ++t0)
      for (int t1 = t0 + 1; t1 < number_of_usable_samples; ++t1)
      {
        double w0 = history_.prev(usable_samples[t0]).w();
        double w1 = history_.prev(usable_samples[t1]).w();
        double spread = utils::square(w0 - w1) + utils::square(w0 - w2_1) + utils::square(w1 - w2_1);
        if (spread > best_spread)
        {
          best_spread = spread;
          i0 = t0;
          i1 = t1;
        }
      }
  }

  Sample const& w1 = history_.prev(usable_samples[i0]);
  Sample const& w0 = history_.prev(usable_samples[i1]);

  double w2_2 = w2_1 * w2_1;
  double w2_3 = w2_2 * w2_1;
  double w2_4 = w2_2 * w2_2;

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
  D <<         w2.dLdw(),
               w1.dLdw(),
               w0.dLdw(),
       w2.Lw() - w1.Lw();

  Eigen::Vector4d C = M.colPivHouseholderQr().solve(D);

  math::Polynomial fourth_degree_approximation(5 COMMA_CWDEBUG_ONLY("w"));
  fourth_degree_approximation[1] = C[0];
  fourth_degree_approximation[2] = C[1];
  fourth_degree_approximation[3] = C[2];
  fourth_degree_approximation[4] = C[3];
  fourth_degree_approximation[0] = w2.Lw() - fourth_degree_approximation(w2_1);
  Dout(dc::notice, "approximation = " << fourth_degree_approximation);

  using namespace cairowindow;

  plot_fourth_degree_approximation_curve_.solve(
      [&fourth_degree_approximation](double w) -> Point { return {w, fourth_degree_approximation(w)}; }, plot_.viewport());
  plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::teal}), plot_fourth_degree_approximation_curve_);

  auto derivative = fourth_degree_approximation.derivative();
  Dout(dc::notice, "derivative = " << derivative);

  plot_derivative_curve_.solve([&derivative](double w) -> Point { return {w, derivative(w)}; }, plot_.viewport());
  plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::magenta}), plot_derivative_curve_);

  double remainder;
  auto quotient = derivative.long_division(w2_1, remainder);
  Dout(dc::notice, "quotient = " << quotient << " with remainder " << remainder);

  std::array<double, 2> zeroes;
  int number_of_zeroes = quotient.get_zeroes(zeroes);

  // It is possible that the zeroes are no usable because they are on the wrong side.
  if (hdirection_ != HorizontalDirection::undecided)
  {
    auto wrong_side = [this, w](double zero) { return (hdirection_ == HorizontalDirection::left) != (zero < w); };
    for (int zero = 0; zero < number_of_zeroes;)
      if (wrong_side(zeroes[zero]))
      {
        if (--number_of_zeroes == 1 && zero == 0)
          zeroes[0] = zeroes[1];
      }
      else
        ++zero;
  }

  if (number_of_zeroes > 1)
    Dout(dc::notice, "with zeroes " << zeroes[0] << " and " << zeroes[1]);
  else if (number_of_zeroes == 1)
    Dout(dc::notice, "with one zero at " << zeroes[0]);
  else
    Dout(dc::notice, "with no zeroes!");

  plot_quotient_curve_.solve([&quotient](double w) -> Point { return {w, 10.0 * quotient(w)}; }, plot_.viewport());
  plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = color::blue}), plot_quotient_curve_);

  if (number_of_zeroes > 0)
  {
    std::array<double, 2> expected_Lw;
    for (int zero = 0; zero < number_of_zeroes; ++zero)
      expected_Lw[zero] = fourth_degree_approximation(zeroes[zero]);
    int best_zero = (number_of_zeroes == 2 &&
        (hdirection_ == HorizontalDirection::right ||
         (hdirection_ == HorizontalDirection::undecided &&
          expected_Lw[1] < expected_Lw[0]))) ? 1 : 0;
    w = zeroes[best_zero];
    expected_Lw_ = expected_Lw[best_zero];
    Dout(dc::notice, history_.current().w() << " --> " << history_.total_number_of_samples() << ": " <<
        w << " [best zero] [expected_Lw: " << expected_Lw_ << "]");
    reset_history();
    Dout(dc::notice(hdirection_ == HorizontalDirection::undecided && number_of_zeroes == 2),
        "Best zero was " << zeroes[best_zero] << " with A(" << zeroes[best_zero] << ") = " <<
        fourth_degree_approximation(zeroes[best_zero]) <<
        " (the other has value " << fourth_degree_approximation(zeroes[1 - best_zero]) << ")");
  }
  else if (hdirection_ == HorizontalDirection::undecided)
  {
    // The "scale" of the (current) parabolic approximation is set to 'edge sample' minus x-coordinate of the vertex (v_x),
    // where 'edge sample' is the sample that is "part of" the parabolic approximation that is the furthest away
    // from the vertex. "Part of" here means that it deviates vertically less than 10% of the vertical distance to
    // the vertex. In that sense we can consider that we "came from" the direction of that edge sample.
    // For example, if the edge sample had x-coordinate w = w_E and we came from the right (w_E > v_x) then
    // scale = w_E - v_x > 0. Subtracting a positive scale thus means we continue in the direction towards to the left.
    // If w_E is on the left of the vertex then scale is negative and subtracting is causes us to continue to the right.
    //
    // Keep going in the same direction.
    w -= new_extreme->approximation().parabola_scale();
    expected_Lw_ = new_extreme->approximation().parabola()(w);
    Dout(dc::notice, history_.current().w() << " --> " << history_.total_number_of_samples() << ": " <<
        w << " [past extreme (no zeroes)] [expected_Lw: " << expected_Lw_ << "]");

    // Note that w was already set to the v_x before, but the test below still works
    // because |v_x - history_.current().w()| was determined to be less than 1% of approximation.parabola_scale().
  }
  else
  {
    // Keep going in the same hdirection.
    w += static_cast<int>(hdirection_) * std::abs(new_extreme->approximation().parabola_scale());
    expected_Lw_ = new_extreme->approximation().parabola()(w);
    Dout(dc::notice, history_.current().w() << " --> " << history_.total_number_of_samples() << ": " <<
        w << " [keep going (no zeroes)] [expected_Lw: " << expected_Lw_ << "]");
  }

  if (hdirection_ == HorizontalDirection::undecided)
  {

    // Now that the decision on which hdirection_ we explore is taken, store that decision.
    hdirection_ = w - history_.current().w() < 0.0 ? HorizontalDirection::left : HorizontalDirection::right;
    Dout(dc::notice, "Initializing horizontal direction to " << hdirection_);
  }

  // Remember in which direction we travelled from this extreme.
  new_extreme->explored(hdirection_);

  // The local extreme was handled.
  state_ = IterationState::done;
  return true;
}

void AcceleratedGradientDescent::update_approximation(bool current_is_replacement)
{
  using namespace gradient_descent;

  math::QuadraticPolynomial old_parabola = plot_approximation_parabola_scale_ptr_->get_parabola();
  ScaleUpdate result = approximation_ptr_->add(&history_.current(), approximation_ptr_ != &current_approximation_, current_is_replacement);

  switch (result)
  {
    case ScaleUpdate::first_sample:
      break;
    case ScaleUpdate::initialized:
      plot_approximation_parabola_scale_ptr_->draw_indicators();
      break;
    case ScaleUpdate::towards_vertex:
      plot_approximation_parabola_scale_ptr_->draw_indicators();
      plot_approximation_parabola_scale_ptr_->draw_old_parabola(old_parabola);
      break;
    case ScaleUpdate::away_from_vertex:
      plot_approximation_parabola_scale_ptr_->draw_indicators();
      break;
    case ScaleUpdate::disconnected:
    {
      reset_history();
      result = approximation_ptr_->add(&history_.current(), false, false);
      break;
    }
  }

  Dout(dc::notice, "approximation = " << *approximation_ptr_ <<
      " (" << utils::print_using(*approximation_ptr_, &Approximation::print_based_on) << ")");

  // Draw the parabolic approximation.
  plot_approximation_curve_.solve(
      [this](double w) -> cairowindow::Point { return {w, approximation_ptr_->parabola()(w)}; }, plot_.viewport());
  plot_.add_bezier_fitter(layer_, curve_line_style_({.line_color = cairowindow::color::red}), plot_approximation_curve_);
}

void AcceleratedGradientDescent::handle_single_sample(Weight& w)
{
  double step;
#ifdef CWDEBUG
  char const* algorithm_str;
#endif

  if (small_step_ == 0.0)       // Not defined yet?
  {
    step = learning_rate_ * -history_.current().dLdw();

    // Did we drop into a (local) minimum as a starting point?!
    if (Scale::almost_zero(w, step))
    {
      if (hdirection_ == HorizontalDirection::undecided)
      {
        // In this case we can't do anything else than just make a step in some random direction.
        step = learning_rate_;
#ifdef CWDEBUG
        algorithm_str = "one sample, derivative is zero, hdirection is unknown";
#endif
      }
      else
      {
        step = static_cast<int>(hdirection_) * learning_rate_;
#ifdef CWDEBUG
        algorithm_str = "one sample, derivative is zero";
#endif
      }
    }
    else if (hdirection_ == HorizontalDirection::undecided)
    {
      // Just gradient descent: move downhill.
#ifdef CWDEBUG
      algorithm_str = "one sample, gradient descent";
#endif
    }
    else
    {
      // Make a step in the same horizontal direction.
      step = static_cast<int>(hdirection_) * std::abs(step);
#ifdef CWDEBUG
      algorithm_str = "one sample, same direction";
#endif
    }
  }
  else
  {
    // small_step_ should only be set once hdirection_ has been set.
    ASSERT(hdirection_ != HorizontalDirection::undecided);
    step = static_cast<int>(hdirection_) * small_step_;
  }

  // This step could still be too small.
  if (approximation_ptr_->parabola_scale().negligible(step))
  {
    // This wouldn't be good because then the new sample will replace
    // the current one and we'd still have just one sample.
    step = approximation_ptr_->parabola_scale().make_significant(step);
#ifdef CWDEBUG
    algorithm_str = "avoiding replacement";
#endif
  }
  w += step;
  expected_Lw_ = approximation_ptr_->parabola()(w);
  Dout(dc::notice, std::setprecision(12) << history_.current().w() << " --> " << history_.total_number_of_samples() <<
      ": " << w << " [" << algorithm_str << "] [expected_Lw:" << expected_Lw_ << "]");

  state_ = IterationState::done;
}

void AcceleratedGradientDescent::handle_parabolic_approximation(Weight& w)
{
  DoutEntering(dc::notice, "AcceleratedGradientDescent::handle_parabolic_approximation(" << w << ")");

  Approximation& approximation(*approximation_ptr_);
  // We have two relevant samples, and thus a parabolic approximation.
  //
  // Case:
  //    A        B        C       D
  // \     /  \     /    ˏ-ˎ     ˏ-ˎ
  // ↑\   /    \   /↑   /   \   /   \
  //   `-´      `-´    /↑    \ /    ↑\
  //
  //                     Jump to the vertex if:
  //      extreme_type   vdirection
  //   A   down             down
  //   B   down             down
  //   C   up               up
  //   D   up               up
  //
  // If the extreme of the current parabolic approximation matches
  // the value of vdirection_ (that the vertical direction that we
  // want to find the next extreme in) then we will jump to that
  // vertex; hdirection_ is not relevant in that case: that is only
  // used when vdirection doesn't match the extreme.

  Sample const& current = approximation.current();
  Sample const& prev = approximation.prev();
  VerticalDirection extreme_type = approximation.has_maximum() ? VerticalDirection::up : VerticalDirection::down;

  if (vdirection_ != extreme_type)
  {
    double v_x = approximation.parabola().vertex_x();
    // If hdirection matches the parabola; jump to the vertex before adding the scale to w.
    if (hdirection_ == (v_x < w ? HorizontalDirection::left : HorizontalDirection::right))
      w = v_x;
    // There is no extreme in the direction that we're going.
    // Just keep going in the same direction as before.
    double step = static_cast<int>(hdirection_) * std::abs(approximation.parabola_scale());
    w += step;
    expected_Lw_ = approximation.parabola()(w);
    Dout(dc::notice, history_.current().w() << " --> " << history_.total_number_of_samples() <<
        ": " << w << " [continue same direction] [expected_Lw:" << expected_Lw_ << "]");
    Dout(dc::notice, (hdirection_ == HorizontalDirection::left ? "Incremented" : "Decremented") << " w with scale (" << std::abs(step) << ")");
    // Abort if the result requires more energy than we have.
    state_ = IterationState::check_energy;
    return;
  }

  // Set w to the value where the derivative of this parabolic approximation is zero.
  Dout(dc::notice, "Setting new_sample to the extreme of parabolic approximation:");
  bool looking_for_maximum = approximation.has_maximum();
  auto ignore = [&current, looking_for_maximum](Sample const& sample) {
    double w0 = current.w();
    double Lw0 = current.Lw();
    double dLdw0 = current.dLdw();
    double w1 = sample.w();
    double Lw1 = sample.Lw();
    double dLdw1 = sample.dLdw();

    // See clamp_check.cxx.
    double dw = w0 - w1;
    double dw3 = std::pow(dw, 3.0);
    double d = (-2.0 * (Lw0 - Lw1) + dw * (dLdw0 + dLdw1)) / dw3;
    double c = (dLdw0 - dLdw1) / (2.0 * dw) - 1.5 * (w0 + w1) * d;
    double b = (w0 * dLdw1 - w1 * dLdw0) / dw + 3.0 * w0 * w1 * d;

    double D = utils::square(c) - 3.0 * b * d;
    if (D >= 0.0)
    {
      double zero = (-c + (looking_for_maximum ? -1.0 : 1.0) * std::sqrt(D)) / (3.0 * d);
      if (std::min(w0, w1) < zero && zero < std::max(w0, w1))
        return false;
    }

    // There is no minimum/maximum in between w0 and w1, therefore
    // do not clamp on this sample: ignore it.
    return true;
  };
  double new_w = history_.clamp(w, approximation.parabola().vertex_x(), ignore, clamped_history_index_);
  if (!clamped_history_index_.undefined())
  {
    state_ = IterationState::clamped;
    return;
  }

  double step = w - new_w;
  w = new_w;
  expected_Lw_ = approximation.parabola().vertex_y();
  state_ = IterationState::vertex_jump;

  double abs_step = std::abs(step);
  Dout(dc::notice, "abs_step = " << abs_step << " (between " <<
      (history_.total_number_of_samples() - 1) << " and " << history_.total_number_of_samples() << " (to be added))");

  // Did we reach the (local) extreme?
  if (abs_step < (extreme_type == VerticalDirection::down ? 0.01 : 0.05) * approximation.parabola_scale())
  {
    Dout(dc::notice, (extreme_type == VerticalDirection::down ? "Minimum" : "Maximum") << " reached: " << abs_step <<
        " < " << (extreme_type == VerticalDirection::down ? 0.01 : 0.05) << " * " << approximation.parabola_scale());
    // If we do not already have at least three relevant samples, then delay reporting a local extreme
    // until after adding one more sample, taking take it won't replace a previous one.
    if (history_.relevant_samples() < 3)
    {
      // Instead of returning this extreme, do a little overshoot.
      double step = static_cast<int>(hdirection_) * Scale::epsilon;
      w += approximation.parabola_scale().make_significant(step);
      expected_Lw_ = approximation.parabola()(w);
      Dout(dc::notice, history_.current().w() << " --> " << history_.total_number_of_samples() << ": " <<
          w << " [jump to vertex + small overshoot] [expected_Lw:" << expected_Lw_ << "]");
      return;
    }
    state_ = IterationState::local_extreme;
  }

  Dout(dc::notice, history_.current().w() << " --> " << history_.total_number_of_samples() << ": " <<
      w << " [jump to vertex] [expected_Lw:" << expected_Lw_ << "]");
}

bool AcceleratedGradientDescent::handle_abort_hdirection(Weight& w)
{
  using namespace gradient_descent;

  // Has a minimum been found at all?
  if (best_minimum_ == extremes_.end())
  {
    Dout(dc::notice, "Aborting going " << hdirection_ << " and terminating search because no extreme has has been found at all!");
    return false;
  }

  // Jump back to the best minimum and continue in the opposite hdirection.
  Dout(dc::notice, "Aborting exploring " << hdirection_ << " of the minimum at " << best_minimum_->vertex_sample() << ".");

  // Was this minimum already explored in both directions?
  if (best_minimum_->done())
    return false;

  ASSERT(hdirection_ != HorizontalDirection::undecided);

  // Change hdirection.
  hdirection_ = opposite(hdirection_);
  Dout(dc::notice, "Changed horizontal direction to " << hdirection_);

  // Restore the current sample and scale to the values belonging to this minimum.
  w = best_minimum_->vertex_sample().w();
  Dout(dc::notice, "Restored w to best minimum at " << w);

  // Do a step with a size equal to the scale in this minimum: we expect it not to drastically change before that point.
  w += static_cast<int>(hdirection_) * std::abs(best_minimum_->approximation().parabola_scale());
  expected_Lw_ = best_minimum_->approximation().parabola()(w);
  Dout(dc::notice, history_.current().w() << " --> " << history_.total_number_of_samples() <<
      ": " << w << " [best minimum, opposite direction] [expected_Lw:" << expected_Lw_ << "]");
  best_minimum_->explored(hdirection_);
  vdirection_ = VerticalDirection::up;
  Dout(dc::notice, "vdirection_ is set to " << vdirection_ << " because we just jumped to a minimum.");

  // Restore the energy to what it was when this minimum was stored.
  energy_.set(best_minimum_->energy(), best_minimum_->vertex_sample().Lw());

  // Restore small_step_ to what is was.
  small_step_ = best_minimum_->plot_approximation_parabola_scale_.scale();

  reset_history();

  // w was successfully updated.
  Dout(dc::notice, "Returning: " << w);
  state_ = IterationState::done;
  return true;
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  Function L;

  double const w_0 = Function::w_0;
  double const w_min = Function::w_min;
  double const w_max = Function::w_max;

  int const steps = 100;

  try
  {
    using namespace cairowindow;
    using namespace gradient_descent;
    using Window = cairowindow::Window;

    // Create a window.
    Window window("Gradient descent of " + L.as_string(), 2*600, 2*450);

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
//      L_min = -100;
//      L_max = 60.0;
    }

    // Create and draw plot area.
    plot::Plot plot(window.geometry(), { .grid = {.color = color::orange} },
        "(1−σ(3(w+55))) (14.6+3.15w+0.451w²) + σ(3(w+55))(0.012 exp((55-w)/10) sin(0.5w+1.92)+1878)", {},
        "w", {},
        "L", {});
    plot.set_xrange({w_min, w_max});
    plot.set_yrange({L_min, L_max});
    plot.add_to(background_layer, false);

    draw::PointStyle point_style({.color_index = 31, .filled_shape = 1});
    draw::LineStyle curve_line_style({.line_width = 1.0});
    draw::TextStyle label_style({.position = draw::centered_left_of, .font_size = 18.0, .offset = 10});
    draw::LineStyle derivative_line_style({.line_color = color::turquoise, .line_width = 1.0});

#if !USE_SLIDERS
    BezierFitter L_fitter([&L](double w) -> Point { return {w, L(w)}; }, plot.viewport());
    auto plot_curve = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(L_fitter));
#endif

    AcceleratedGradientDescent agd(0.1, L_max, plot, second_layer, point_style, label_style);

#if USE_SLIDERS
    // amplitude = 0.012027, level = 1878.38, phase = 1.91892
    draw::TextStyle slider_style({.position = draw::centered_below, .font_size = 18.0, .offset = 10});
    auto slider_amplitude = plot.create_slider(second_layer, {978, 83, 7, 400}, 0.012027, 0.01, 0.04);
    auto slider_amplitude_label = plot.create_text(second_layer, slider_style, Pixel{978, 483}, "amplitude");
    auto slider_level = plot.create_slider(second_layer, {1028, 83, 7, 400}, 1878.38, 1500, 2500);
    auto slider_level_label = plot.create_text(second_layer, slider_style, Pixel{1028, 483}, "level");
    auto slider_phase = plot.create_slider(second_layer, {1078, 83, 7, 400}, 1.91892, 0.0, 2.0 * M_PI);
    auto slider_phase_label = plot.create_text(second_layer, slider_style, Pixel{1078, 483}, "phase");
#endif

    // Loop over iterations of w.
    for (Weight w(w_0);;)
    {
      // Suppress immediate updating of the window for each created item, in order to avoid flickering.
      window.set_send_expose_events(false);

#if USE_SLIDERS
      L.set_sliders(slider_amplitude.value(), slider_level.value(), slider_phase.value());
      Dout(dc::notice,
          "amplitude = " << slider_amplitude.value() << ", level = " << slider_level.value() << ", phase = " << slider_phase.value());

      BezierFitter L_fitter([&L](double w) -> Point { return {w, L(w)}; }, plot.viewport());
      auto plot_curve = plot.create_bezier_fitter(second_layer, curve_line_style, std::move(L_fitter));
#else
      // Replace w with a new point until the global minimum has been reached.
      if (!agd(w, L(w), L.derivative(w)))
        break;
#endif

      // Flush all expose events related to the drawing done above.
      window.set_send_expose_events(true);

      // Block until a key is pressed.
      if (!window.handle_input_events())
        break;          // Program must be terminated.

      Dout(dc::notice, "------------------------------------");
    }

#if !USE_SLIDERS
    if (agd.success())
      Dout(dc::notice, "Found global minimum " << agd.minimum());
#endif

    event_loop.join();
  }
  catch (AIAlert::Error const& error)
  {
    Dout(dc::warning, error);
  }

  Dout(dc::notice, "Leaving main()");
}
