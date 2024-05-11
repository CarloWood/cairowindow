#include "sys.h"
#include "Polynomial.h"
#include "QuadraticPolynomial.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Plot.h"
#include "cairowindow/symbolic/symbolic.h"
#include "utils/square.h"
#include "utils/AIAlert.h"
#include "utils/debug_ostream_operators.h"
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

#if 1
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
#else
constexpr double a2 = 14.6;
constexpr double b2 = 3.15;
constexpr double c2 = 0.451;

class Function
{
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

class Sample
{
 private:
  double w_;                    // w
  mutable bool initialized_;    // Set to true iff Lw_ and dLdw_ correspond to the current w_ (delayed initialization).
  mutable double Lw_;           // Cache of L(w).
  mutable double dLdw_;         // Cache of L'(w).

  cairowindow::plot::Point P_;
  cairowindow::plot::Text P_label_;

 public:
  // Create an uninitialized Sample.
  Sample() = default;

  // Create Sample for which only w is known.
  Sample(double w) : w_(w), initialized_(false) { }

  Sample(double w, double Lw, double dLdw, cairowindow::plot::Point const& point, cairowindow::plot::Text const& label) :
    w_(w), Lw_(Lw), dLdw_(dLdw), initialized_(true), P_(point), P_label_(label) { }

  Sample& operator-=(double delta_w)
  {
    w_ -= delta_w;
    initialized_ = false;
    return *this;
  }

  Sample& operator+=(double delta_w)
  {
    w_ += delta_w;
    initialized_ = false;
    return *this;
  }

  std::string const& label() const
  {
    return P_label_.text();
  }

  double w() const
  {
    return w_;
  }

  double Lw() const
  {
    ASSERT(initialized_);
    return Lw_;
  }

  double dLdw() const
  {
    ASSERT(initialized_);
    return dLdw_;
  }

  double Lw(Function const& L) const
  {
    if (!initialized_)
    {
      Lw_ = L(w_);
      dLdw_ = L.derivative(w_);
      initialized_ = true;
    }
    return Lw_;
  }

  double dLdw(Function const& L) const
  {
    if (!initialized_)
    {
      Lw_ = L(w_);
      dLdw_ = L.derivative(w_);
      initialized_ = true;
    }
    return dLdw_;
  }

  void set_values(double w, double Lw, double dLdw)
  {
    w_ = w;
    Lw_ = Lw;
    dLdw_ = dLdw;
  }

  void set_point(cairowindow::plot::Point const& P)
  {
    P_ = P;
  }

  void set_label(cairowindow::plot::Text const& text)
  {
    P_label_ = text;
  }

  cairowindow::Point const& P() const { return P_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << label() << " (at " << w_ << ")";
  }
#endif
};

class Scale
{
 public:
  static constexpr double epsilon = 1e-30;
  static cairowindow::draw::ConnectorStyle const s_indicator_style;

 private:
  double scale_{};                      // An indication of what changes to w are significant.
  bool has_sample_{};                   // True iff the below is value.
  double edge_sample_w_;                // The value of w that corresponds to this scale: edge_sample_w_ - scale_
                                        // should be more or less equal to the vertex of the parabola_.
  double edge_sample_Lw_;               // Cached value of L(edge_sample_w_). Should also be more or less equal
                                        // to the value of the parabola_ at edge_sample_w_.
  math::QuadraticPolynomial parabola_;  // The last (previous) second degree polynomial fit (passed to initialize/update).

  // Used to visualize the Scale:
  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> layer_;
  cairowindow::plot::Connector plot_indicator_;
  cairowindow::plot::Text plot_scale_text_;
  cairowindow::plot::Line plot_vertical_line_through_w_;
  cairowindow::plot::Line plot_vertical_line_through_v_;

  // Temporary curves (used while developing this class).
  cairowindow::plot::BezierFitter plot_old_parabola_;

 private:
  void draw_indicators()
  {
    ASSERT(has_sample_);

    double x1 = parabola_.vertex_x();
    double x2 = edge_sample_w_;
    double scale_y = plot_.yrange().min() + 0.5 * plot_.yrange().size();
    plot_indicator_ = cairowindow::plot::Connector{{x1, scale_y}, {x2, scale_y},
        cairowindow::Connector::open_arrow, cairowindow::Connector::open_arrow};
    plot_.add_connector(layer_, s_indicator_style, plot_indicator_);
    plot_scale_text_ = plot_.create_text(layer_, {{.position = cairowindow::draw::centered_above, .offset = 2.0}},
        cairowindow::Point{(x1 + x2) / 2, scale_y}, "scale");

    plot_vertical_line_through_v_ = cairowindow::plot::Line{{parabola_.vertex_x(), scale_y}, cairowindow::Direction::up};
    plot_.add_line(layer_, s_indicator_style({.dashes = {3.0, 3.0}}), plot_vertical_line_through_v_);

    plot_vertical_line_through_w_ = cairowindow::plot::Line{{edge_sample_w_, scale_y}, cairowindow::Direction::up};
    plot_.add_line(layer_, s_indicator_style, plot_vertical_line_through_w_);
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

 public:
  Scale(cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer) : plot_(plot), layer_(layer)
  {
  }

  Scale& operator=(Scale const& scale)
  {
    scale_ = scale.scale_;
    has_sample_ = scale.has_sample_;
    if (has_sample_)
    {
      edge_sample_w_ = scale.edge_sample_w_;
      edge_sample_Lw_ = scale.edge_sample_Lw_;
      parabola_ = scale.parabola_;
    }

    erase_indicators();
    return *this;
  }

  void reset(double scale)
  {
    scale_ = scale;
    has_sample_ = false;
    erase_indicators();
  }

  operator double() const { ASSERT(scale_ != 0.0); return std::abs(scale_); }
  double or_zero() const { /* Never return really zero */ return std::max(epsilon, scale_); }

  static bool almost_zero(double w, double step)
  {
    return std::abs(step) < epsilon || std::abs(step) < 1e-6 * std::abs(w);
  }

  void initialize(Sample const& prev, Sample const& current, math::QuadraticPolynomial const& parabola)
  {
    DoutEntering(dc::notice, "Scale::initialize(" << prev << ", " << current << ", " << parabola << ")");
    // Only call initialize once.
    ASSERT(!has_sample_);
    double const v = parabola.vertex_x();
    Dout(dc::notice, "v = " << v);
    // Not sure it can happen that current is further away, but in case it does do this test.
    // We want to set the scale_ to the largest value that still makes sense: the distance from
    // the sample (that participated in creating this parabolic fit, that is, prev and current)
    // that is the furthest away from the vertex.
    Sample const& edge_sample = (std::abs(prev.w() - v) > std::abs(current.w() - v)) ? prev : current;
    scale_ = edge_sample.w() - v;
    Dout(dc::notice, "scale was set to " << scale_);
    edge_sample_w_ = edge_sample.w();
    edge_sample_Lw_ = edge_sample.Lw();
    parabola_ = parabola;
    has_sample_ = true;

    draw_indicators();
  }

  void update(Sample const& prev, Sample const& current, math::QuadraticPolynomial const& parabola)
  {
    DoutEntering(dc::notice, "Scale::update(" << prev << ", " << current << ", " << parabola << ")");

    // Call initialize if there is nothing to update because the Scale was reset.
    if (!has_sample_)
    {
      initialize(prev, current, parabola);
      return;
    }

    // Get the x coordinate (w value) of the vertex of the new parabola.
    double const v_x = parabola.vertex_x();
    // Pick the sample (from prev and current) that is horizontally the furthest away from the vertex.
    Sample const& edge_sample = (std::abs(prev.w() - v_x) > std::abs(current.w() - v_x)) ? prev : current;
    // Get the y-coordinate at the w value of the stored edge sample, according to the new parabola.
    double const new_Lw_stored_edge_sample = parabola(edge_sample_w_);
    // Get the vertical distance from the stored edge sample to the new parabola.
    double const distance_stored_edge_sample_to_parabola = std::abs(new_Lw_stored_edge_sample - edge_sample_Lw_);
    // Get the vertical distance from the stored edge sample to the vertex of the new parabola (equal to abs(parabola.height(edge_sample_w_))).
    double const abs_stored_edge_sample_height = std::abs(new_Lw_stored_edge_sample - parabola.vertex_y());
    // If the stored sample vertically deviates more than 10%, we need to adjust the stored sample values.
    if (distance_stored_edge_sample_to_parabola > 0.1 * abs_stored_edge_sample_height)
    {
      Dout(dc::notice, "Discarded old sample at L(" << edge_sample_w_ << ") = " << edge_sample_Lw_ << " because " <<
          distance_stored_edge_sample_to_parabola << " > 0.1 * " << abs_stored_edge_sample_height);

      static constexpr int left = -1;
      static constexpr int right = 1;
      // On which side of the vertex is the old sample?
      int hside = edge_sample_w_ < v_x ? left : right;

      int count;
      std::array<double, 4> toggles;
      bool close = parabola.equal_intervals(parabola_, toggles, count);
      double new_edge_w;

      // Run over intervals from left-to-right when hside == left, and from right-to-left when hside == right.
      //
      //          0     1     2     3
      //    -inf  |     |     |     |  +inf     <-- interval
      //   begin <0    <1    <2    <3  end      <-- running from left-to-right
      //     end  0>    1>    2>    3> begin    <-- running from right-to-left
      //
      //                   ^     ^        ^
      //                   |     |        |
      //                  v_x    |      edge_sample_w_
      //                    edge_sample.w()
      //
      if (hside == left)
      {
        // Find the interval that edge_sample_w_ is in.
        int i = 0;
        while (i != count)
        {
          if (edge_sample_w_ < toggles[i])
            break;
          ++i;
          close = !close;
        }
        // Use the old value, unless it is too far away from the new parabola. Don't go beyond the new edge sample.
        new_edge_w = std::min(close ? edge_sample_w_ : toggles[i], edge_sample.w());
      }
      else
      {
        // Find the interval that edge_sample_w_ is in.
        int i = count - 1;
        while (i != -1)
        {
          if (edge_sample_w_ > toggles[i])
            break;
          --i;
          close = !close;
        }
        // Use the old value, unless it is too far away from the new parabola. Don't go beyond the new edge sample.
        new_edge_w = std::max(close ? edge_sample_w_ : toggles[i], edge_sample.w());
      }

      edge_sample_w_ = new_edge_w;
      edge_sample_Lw_ = utils::almost_equal(new_edge_w, edge_sample.w(), 1e-3) ? edge_sample.Lw() : parabola(new_edge_w);
    }
    scale_ = edge_sample_w_ - v_x;
    Dout(dc::notice, "scale was set to " << scale_);

    // Draw the old parabola, for debugging purposes.
    using namespace cairowindow;
    plot_old_parabola_.solve([&](double w) -> Point { return {w, parabola_(w)}; }, plot_.viewport());
    plot_.add_bezier_fitter(layer_, {{.line_color = color::light_red, .line_width = 1.0}}, plot_old_parabola_);

    parabola_ = parabola;

    draw_indicators();
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << scale_;
  }
#endif
};

class SampleMinimum : public Sample
{
 private:
  Scale scale_;
  double energy_;

 public:
  SampleMinimum(Sample const& minimum, Scale const& scale, double energy) : Sample(minimum), scale_(scale), energy_(energy) { }

  Scale const& scale() const { return scale_; }
  double energy() const { return energy_; }
};

class History
{
 public:
  static constexpr int size = 9;

 private:
  Function const& L_;
  std::array<Sample, size> samples_;
  int prev_ = -1;
  int old_samples_ = 0;                 // Samples elsewhere that should not be used for fitting a polynomial.
  int total_number_of_samples_ = 0;     // The number of samples taken (not necessarily equal to the number that is (still) in the history).
  int current_ = -1;

  cairowindow::plot::Plot& plot_;
  boost::intrusive_ptr<cairowindow::Layer> const& layer_;
  cairowindow::draw::PointStyle const& point_style_;
  cairowindow::draw::TextStyle const& label_style_;

 public:
  History(Function const& L, cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer,
      cairowindow::draw::PointStyle const& point_style, cairowindow::draw::TextStyle const& label_style) :
    L_(L), plot_(plot), layer_(layer), point_style_(point_style), label_style_(label_style) { }

  void add(Sample& sample, Scale const& scale)
  {
    double w = sample.w();
    double Lw = sample.Lw(L_);
    double dLdw = sample.dLdw(L_);
    // If the new sample is very close to the current one, don't add it, just replace the current sample.
    bool almost_equal = (current_ != -1 && std::abs(samples_[current_].w() - w) < 0.001 * scale.or_zero());
    if (almost_equal)
    {
      Dout(dc::notice, "Replacing history point " << (total_number_of_samples_ - 1) << " at w = " << samples_[current_].w() << " with " << w);
      samples_[current_].set_values(w, Lw, dLdw);
      samples_[current_].set_point(plot_.create_point(layer_, point_style_, {w, Lw}));
    }
    else
    {
      prev_ = current_;
      current_ = (prev_ + 1) % size;
      sample.set_point(plot_.create_point(layer_, point_style_, {w, Lw}));
      sample.set_label(plot_.create_text(layer_, label_style_({.position = cairowindow::draw::centered_below}), sample.P(),
            std::to_string(total_number_of_samples_)));
      Dout(dc::notice, "Appended to history: point " << total_number_of_samples_ << " at w = " << sample.w());
      samples_[current_] = sample;
      ++total_number_of_samples_;
    }
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
    if (closest_index != current_)
    {
      prev_ = current_;
      current_ = (prev_ + 1) % size;
      if (current_ != closest_index)
        samples_[current_] = samples_[closest_index];
      samples_[current_].set_label(plot_.create_text(
            layer_, label_style_({.position = cairowindow::draw::centered_above}), samples_[closest_index].P(),
            std::to_string(total_number_of_samples_)));
      Dout(dc::notice, "Created duplicate of " << closest_index <<
          " with label " << total_number_of_samples_ << " at w = " << samples_[current_].w());
      old_samples_ = total_number_of_samples_;
      ++total_number_of_samples_;
    }
    else
    {
      Dout(dc::notice, "Current sample " << (total_number_of_samples_ - 1) << " at w = " << samples_[current_].w() << " is already closest.");
      old_samples_ = total_number_of_samples_ - 1;
    }
  }

  int relevant_samples() const
  {
    return total_number_of_samples_ - old_samples_;
  }

  static int before(int i) { ASSERT(0 <= i); return (i + size - 1) % size; }

  Sample const& current() const { ASSERT(current_ != -1); ASSERT(current_ < total_number_of_samples_); return samples_[current_]; }
  Sample const& prev() const { ASSERT(prev_ != -1); ASSERT(prev_ < total_number_of_samples_); return samples_[prev_]; }
  Sample const& prev_prev() const { ASSERT(prev_ != -1); ASSERT(prev_ < total_number_of_samples_); return samples_[before(prev_)]; }

  int total_number_of_samples() const { return total_number_of_samples_; }
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

  double const w_min = -20.0;
  double const w_max = 80.0;
  int const steps = 100;

  double const w_0 = -5.0;
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
//      L_max = 1.25;
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
    History history(L, plot, second_layer, point_style, label_style);

    // Initial value.
    Sample new_sample(w_0);
    Dout(dc::notice, "Initial value of new_sample; w = " << new_sample.w());
    double current_learning_rate = learning_rate;
    Scale scale(plot, second_layer);
    Dout(dc::notice, "Initial value of scale is " << scale);
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

    while (true)
    {
      // Add new sample to the history.
      history.add(new_sample, scale);
      Sample const& current = history.current();        // Currently, the same as new_sample, but more clear that it
                                                        // was already added to the history. If new_sample is updated
                                                        // then current is not, of course.

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

      enum class StepKind
      {
        keep_going,             // Just subtract current_learning_rate times derivative.
        done,                   // new_sample was already updated.
        abort                   // Stop going in the current vdirection.
      };

      // The default is that new_sample is set to the extreme of some polynomial.
      StepKind step_kind = StepKind::done;

      //===================================================================================================
      // Create a parabolic approximation from the last two samples (or a line if we only have one sample).
      // If we have only one sample, then new_sample and energy accordingly.

      math::QuadraticPolynomial parabolic_approximation;

      // Do we have only one sample?
      if (history.total_number_of_samples() == 1)
      {
        double const dLdw = current.dLdw(L);
        double const w = current.w();

        // If we have just one point, then the approximation is a linear function:
        //
        // A(w) = coef[0] + L'(w) w
        parabolic_approximation[1] = dLdw;

        // Did we drop into a (local) minimum as a starting point?!
        if (Scale::almost_zero(w, dLdw))
        {
          // In this case we can't do anything else than just make a step in some random vdirection.
          new_sample -= learning_rate;
        }
        else
        {
          // Just gradient descent: move downhill.
          new_sample -= learning_rate * dLdw;
        }
        // Force an update.
        energy.update(new_sample.Lw(L));
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
          inverse_det * (2.0 * current.Lw(L) - 2.0 * prev.Lw(L) + (prev.dLdw() - current.dLdw()) * (current.w() + prev.w()));
      }
      parabolic_approximation[0] = current.Lw(L) - parabolic_approximation(current.w());

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
          // Just keep going in the same vdirection as before, but accelerating
          // (increase the learning rate with a factor of 1.4).
          current_learning_rate *= 1.4;
          step_kind = StepKind::keep_going;
        }
        else
        {
          // Set w to the value where the derivative of this parabolic approximation is zero.
          Dout(dc::notice, "Setting new_sample to the extreme of parabolic approximation:");
          double step = beta_inverse * current.dLdw();
          new_sample -= step;
          Dout(dc::notice, history.current().w() << " --> " << history.total_number_of_samples() << ": " << new_sample.w());

          // Initialize or update the scale with the new parabolic approximation.
          if (relevant_samples == 2)
            scale.initialize(prev, current, parabolic_approximation);
          else
            scale.update(prev, current, parabolic_approximation);
          // Update the current kinetic energy. If this is an overshoot, just increase the energy to 0.
          energy.update(new_sample.Lw(L));

          // With the new extreme insight, adjust the learning rate to the new scale.
          current_learning_rate = 0.1 * std::abs(beta_inverse);

          if (relevant_samples > 2)
          {
            double abs_step = std::abs(step);
            Dout(dc::notice, "abs_step = " << abs_step << " (between " <<
                (history.total_number_of_samples() - 1) << " and " << history.total_number_of_samples() << " (to be added))");

            // Did we reach the (local) extreme?
            if (abs_step < 0.01 * scale)
            {
              double w2_1 = new_sample.w();
              double w2_2 = w2_1 * w2_1;
              double w2_3 = w2_2 * w2_1;
              double w2_4 = w2_2 * w2_2;

              // If the new sample is too close to the current one, then ignore current.
              bool skip_sample = std::abs(w2_1 - current.w()) < 0.001 * scale;
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
              D <<                new_sample.dLdw(L),
                                           w1.dLdw(),
                                           w0.dLdw(),
                          new_sample.Lw(L) - w1.Lw();

              Eigen::Vector4d C = M.colPivHouseholderQr().solve(D);

              math::Polynomial approximation(5 COMMA_CWDEBUG_ONLY("w"));
              approximation[1] = C[0];
              approximation[2] = C[1];
              approximation[3] = C[2];
              approximation[4] = C[3];
              approximation[0] = new_sample.Lw() - approximation(new_sample.w());
              Dout(dc::notice, "approximation = " << approximation);

              plot_approximation_curve.solve([&approximation](double w) -> Point { return {w, approximation(w)}; }, plot.viewport());
              plot.add_bezier_fitter(second_layer, curve_line_style({.line_color = color::teal}), plot_approximation_curve);

              auto derivative = approximation.derivative();
              Dout(dc::notice, "derivative = " << derivative);

              plot_derivative_curve.solve([&derivative](double w) -> Point { return {w, derivative(w)}; }, plot.viewport());
              plot.add_bezier_fitter(second_layer, curve_line_style({.line_color = color::magenta}), plot_derivative_curve);

              double remainder;
              auto quotient = derivative.long_division(new_sample.w(), remainder);
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

              scale.reset(number_of_zeroes < 2 ? scale : zeroes[1] - zeroes[0]);

              // Did we reach the (local) extreme?
              if (abs_step < 0.01 * scale)
              {
                Dout(dc::notice, (vdirection == up ? "Maximum" : "Minimum") << " reached: " << abs_step << " < 0.01 * " << scale);

                // Store the found extreme in the history!
                history.add(new_sample, scale);

                // Change vdirection.
                vdirection = -vdirection;
                Dout(dc::notice, "vdirection is set to " << (vdirection == up ? "up" : "down") <<
                    "; hdirection = " << (hdirection == left ? "left" : "right"));

                if (number_of_zeroes > 0)
                {
                  int best_zero = (number_of_zeroes == 2 &&
                      (hdirection == right || (hdirection == unknown && approximation(zeroes[1]) < approximation(zeroes[0])))) ? 1 : 0;
                  new_sample = zeroes[best_zero];
                  Dout(dc::notice, "Set w to found extreme: w = " << new_sample.w());
                  Dout(dc::notice(hdirection == unknown && number_of_zeroes == 2),
                      "Best zero was " << zeroes[best_zero] << " with A(" << zeroes[best_zero] << ") = " <<
                      approximation(zeroes[best_zero]) << " (the other has value " << approximation(zeroes[1 - best_zero]) << ")");
                  // Update the current kinetic energy. If this is an overshoot, just increase the energy to 0.
                  energy.update(new_sample.Lw(L));
                }

                if (hdirection == unknown)
                {
                  if (number_of_zeroes == 0)
                  {
                    // Keep going in the same vdirection.
                    new_sample -= scale;
                    energy.update(new_sample.Lw(L));
                  }

                  // Now that the decision on which hdirection we explore is taken, store that decision.
                  hdirection = new_sample.w() - history.current().w() < 0.0 ? left : right;
                  Dout(dc::notice, "Set horizontal direction to " << (hdirection == left ? "left" : "right"));
                }

                // Store this extreme.
                std::list<SampleMinimum>::iterator new_extreme;
                if (hdirection == right)
                  new_extreme = extremes.emplace(extremes.end(), history.current(), scale, energy.energy());
                else
                  new_extreme = extremes.emplace(extremes.begin(), history.current(), scale, energy.energy());

                if (vdirection == up)
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

                if (step_kind == StepKind::done)
                {
                  // Re-add old sample again to the history; add one closest to the target.
                  // This also resets the number of relevant samples to 1.
                  history.append_closest_to(new_sample.w());
                }
              }
            }
          }
        }
      }

      if (step_kind == StepKind::keep_going)
      {
        // There is no new extreme insight yet; just keep going (up or) down hill.
        //FIXME: multiplying with the derivative makes little sense in the "keep_going" case!
        double step = current_learning_rate * history.current().dLdw();
        double w = new_sample.w();
        if (vdirection == down)
          new_sample -= step;
        else
          new_sample += step;
        bool success = energy.maybe_update(new_sample.Lw(L));
        if (!success)
        {
          new_sample = w;
          step_kind = StepKind::abort;
        }
        else
          Dout(dc::notice, ((vdirection == (step > 0.0 ? up : down)) ? "Incremented" : "Decremented") <<
              " new_sample with learning rate of " << current_learning_rate << " and slope " <<
              history.current().dLdw() << " with " << std::abs(step));
      }

      if (step_kind == StepKind::abort)
      {
        // Jump back to the best minimum and continue in the opposite hdirection.

        // Has a minimum been found at all?
        if (best_minimum == extremes.end())
          break;

        // Restore the current sample and scale to the values belonging to this minimum.
        new_sample = *best_minimum;
        scale = best_minimum->scale();

        // Forget any acceleration that happened before.
        current_learning_rate = learning_rate;

        // Change hdirection.
        hdirection = -hdirection;
        Dout(dc::notice, "Changed horizontal direction to " << (hdirection == left ? "left" : "right"));

        // Do a step with a size equal to the scale in this minimum: we expect it to be a parabola up till that point.
        new_sample += hdirection * std::abs(scale);

        energy.set(best_minimum->energy(), new_sample.Lw(L));
      }

      if (step_kind != StepKind::done)
        Dout(dc::notice, history.current().w() << " --> " << history.total_number_of_samples() << ": " << new_sample.w());
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
