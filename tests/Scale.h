#pragma once

#include "QuadraticPolynomial.h"
#include "utils/AIRefCount.h"
#include <type_traits>
#include <utility>
#include "debug.h"

namespace gradient_descent {
using utils::has_print_on::operator<<;

class Scale;

template<typename T>
concept ConceptScale = std::is_base_of_v<Scale, T>;

class Scale : public AIRefCount
{
 public:
  static constexpr double epsilon = 1e-30;

 protected:
  double scale_{};                      // An indication of what changes to w are significant.
  bool has_sample_{};                   // True iff the below is value.
  double edge_sample_w_;                // The value of w that corresponds to this scale: edge_sample_w_ - scale_
                                        // should be more or less equal to the vertex of the parabola_.
  double edge_sample_Lw_;               // Cached value of L(edge_sample_w_). Should also be more or less equal
                                        // to the value of the parabola_ at edge_sample_w_.
  math::QuadraticPolynomial parabola_;  // The last (previous) second degree polynomial fit (passed to initialize/update).

 protected:
  // Use create<>.
  Scale() = default;
  Scale(Scale const& scale)  = delete;

 public:
  Scale& operator=(Scale const& scale) = delete;
#if 0
  // Call draw_indicators() after calling this!
  Scale& operator=(Scale const& scale)
  {
    DoutEntering(dc::notice, "Scale::operator=(" << scale << ")");

    scale_ = scale.scale_;
    has_sample_ = scale.has_sample_;
    if (has_sample_)
    {
      edge_sample_w_ = scale.edge_sample_w_;
      edge_sample_Lw_ = scale.edge_sample_Lw_;
      parabola_ = scale.parabola_;
    }
    return *this;
  }
#endif

  template<ConceptScale T, typename... Args>
  static boost::intrusive_ptr<T> create(Args&&... args)
  {
    return new T(std::forward<Args>(args)...);
  }

  void reset(double scale)
  {
    DoutEntering(dc::notice, "Scale::reset()");
    scale_ = scale;
    has_sample_ = false;
  }

  double learning_rate() const
  {
    // This only makes sense when parabola_ is defined.
    ASSERT(has_sample_);

    // Return learning_rate = 0.1 / |beta| where beta = (f'(w1) - f'(w0)) / (w1 - w0), for any two distinct values of w0 and w1.
    // Let the parabola_ be f(x) = a + b x + c x^2 with derivative: f'(x) = b + 2c x.
    // f'(w1) - f'(w0) = b + 2c w1 - (b + 2c w0) = 2(w1 - w0) c.
    // beta = 2(w1 - w0) c / (w1 - w0) = 2c (the constant value of the second derivative).

    // Return 0.1 / beta.
    return 0.05 / std::abs(parabola_[2]);
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
    double const v_x = parabola.vertex_x();
    Dout(dc::notice, "v_x = " << v_x);
    // Not sure it can happen that current is further away, but in case it does do this test.
    // We want to set the scale_ to the largest value that still makes sense: the distance from
    // the sample (that participated in creating this parabolic fit, that is, prev and current)
    // that is the furthest away from the vertex.
    Sample const& edge_sample = (std::abs(prev.w() - v_x) > std::abs(current.w() - v_x)) ? prev : current;
    scale_ = edge_sample.w() - v_x;
    Dout(dc::notice, "scale was set to " << scale_);
    edge_sample_w_ = edge_sample.w();
    edge_sample_Lw_ = edge_sample.Lw();
    parabola_ = parabola;
    has_sample_ = true;
  }

  bool update(Sample const& prev, Sample const& current, math::QuadraticPolynomial const& parabola)
  {
    DoutEntering(dc::notice, "Scale::update(" << prev << ", " << current << ", " << parabola << ")");

    // Call initialize if there is nothing to update because the Scale was reset.
    if (!has_sample_)
    {
      initialize(prev, current, parabola);
      return false;
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

    parabola_ = parabola;

    return true;
  }

  bool update(Sample const& new_sample)
  {
    // If there is no parabola, then there is nothing to update.
    if (!has_sample_)
      return false;

    // Get the x coordinate (w value) of the vertex of the parabola.
    double const v_x = parabola_.vertex_x();
    double new_scale = new_sample.w() - v_x;
    if (std::abs(new_scale) > std::abs(scale_))
    {
      // Get the y-coordinate at the w value of the new sample, according to the parabola.
      double const parabola_Lw_at_new_sample = parabola_(new_sample.w());
      // Get the vertical distance from the new sample to the new parabola.
      double const distance_new_sample_to_parabola = std::abs(parabola_Lw_at_new_sample - new_sample.Lw());
      // Get the vertical distance from the new sample to the vertex of the parabola.
      double const abs_new_sample_height = std::abs(new_sample.Lw() - parabola_.vertex_y());
      // If the new sample vertically deviates less than 10%, we can use the new sample as new edge sample and adjust the scale accordingly.
      if (distance_new_sample_to_parabola < 0.1 * abs_new_sample_height)
      {
        scale_ = new_scale;
        edge_sample_w_ = new_sample.w();
        edge_sample_Lw_ = new_sample.Lw();
        Dout(dc::notice, "scale was set to " << scale_);
        return true;
      }
    }

    return false;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{scale:" << scale_ << ", has_sample:" << std::boolalpha << has_sample_ << "}";
  }
#endif
};

} // namespace gradient_descent
