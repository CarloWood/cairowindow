#pragma once

#include "Sample.h"
#include "HorizontalDirection.h"
#include "../CubicPolynomial.h"
#include "utils/debug_ostream_operators.h"
#include "utils/almost_equal.h"
#include "utils/macros.h"
#include <type_traits>
#include <utility>
#include "debug.h"

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Scale;

template<typename T>
concept ConceptScale = std::is_base_of_v<Scale, T>;

enum class ScaleUpdate
{
  first_sample,                 // When there is only one sample.
  initialized,                  // When there are two samples for the first time.
  towards_vertex,               // When we already had two samples (and therefore an "old" parabolic approximation)
                                // and the new sample is at the vertex of the new parabola.
  away_from_vertex,             // When we already had two samples (and therefore an "old" parabolic approximation)
                                // and the new sample was moved away from the vertex.
  disconnected                  // Same as away_from_vertex but the scale wasn't increased.
};

#ifdef CWDEBUG
inline std::string to_string(ScaleUpdate scale_update)
{
  switch (scale_update)
  {
    AI_CASE_RETURN(ScaleUpdate::first_sample);
    AI_CASE_RETURN(ScaleUpdate::initialized);
    AI_CASE_RETURN(ScaleUpdate::towards_vertex);
    AI_CASE_RETURN(ScaleUpdate::away_from_vertex);
    AI_CASE_RETURN(ScaleUpdate::disconnected);
  }
  AI_NEVER_REACHED
}

inline std::ostream& operator<<(std::ostream& os, ScaleUpdate scale_update)
{
  return os << to_string(scale_update);
}
#endif

enum class CriticalPointType
{
  minimum,
  maximum,
  inflection_point
};

class Scale
{
 public:
  static constexpr double epsilon = 1e-30;
  static constexpr int left = 0;        // Index into the arrays.
  static constexpr int right = 1;

  struct Members {
    bool valid;                 // True iff the members below are valid.
    double scale;               // An indication of what changes to w are significant, relative to the critical point.
                                // This value is always larger than zero.
    double edge_sample_w;       // The value of w that corresponds to scale: edge_sample_w +/- scale_ should be more
                                // or less equal to the critical point (+ if this is left, - if this is right).
    double edge_sample_Lw;      // Cached value of L(edge_sample_w).
                                // Should also be more or less equal to the value of the cubic_(edge_sample_w).
  };

 protected:
  math::CubicPolynomial cubic_;         // The last (previous) third degree polynomial fit (passed to initialize/update).
  CriticalPointType type_;              // Whether distances are measured relative to the minimum, maximum or
                                        // the inflection point (if the cubic has no extremes).
  double critical_point_w_;             // Cached value of the x-coordinate of the critical point.
  std::array<Members, 2> members_{};    // Separate members for values related to samples left and right of the critical point.

 public:
  Scale() = default;

 public:
  Scale& operator=(Scale const& scale) = delete;

  // Accessors.
  math::CubicPolynomial const& cubic() const { return cubic_; }
  CriticalPointType type() const { return type_; }
  bool has_sample(int side) const { return members_[side].valid; }
  double value(int side) const { ASSERT(members_[side].valid); return members_[side].scale; }
  double critical_point_w() const { return critical_point_w_; }
  double edge_sample_w(int side) const { ASSERT(members_[side].valid); return members_[side].edge_sample_w; }
  double edge_sample_Lw(int side) const { ASSERT(members_[side].valid); return members_[side].edge_sample_Lw; }

  double step(double w, HorizontalDirectionToInt direction) const
  {
    int side = (w > critical_point_w_) ? right : left;
    ASSERT(members_[side].valid);
    return direction * members_[side].scale;
  }

  double avg() const
  {
    // The number of valid sides.
    int number_of_valid_sides = (members_[left].valid? 1 : 0) + (members_[right].valid ? 1 : 0);
    // Do not call this function if this Scale has no samples at all yet.
    ASSERT(number_of_valid_sides > 0);
    return ((members_[left].valid ? members_[left].scale : 0.0) + (members_[right].valid ? members_[right].scale : 0.0)) / number_of_valid_sides;
  }

  void reset()
  {
    DoutEntering(dc::notice, "Scale::reset()");
    members_[left].valid = false;
    members_[right].valid = false;
  }

  double learning_rate(int side) const
  {
    // This only makes sense when cubic_ is defined.
    ASSERT(has_sample(side));

    // Return a learning_rate such that, for the given side,
    // * if this is a minimum,
    //     minimum = edge_sample_w - learning_rate * cubic_.derivative(edge_sample_w)
    //     [this is the "normal" gradient descent algorithm with the idea
    //     that if we are left/right of the minimum, the derivative will be negative/positive;
    //     therefore we need to subtract (a positive factor times) the derivative.]
    // * if this is a maximum,
    //     maximum = edge_sample_w + learning_rate * cubic_.derivative(edge_sample_w)
    //     [likewise, if we are left/right of a maximum, the derivative will be positive/negative;
    //     therefore we need to add (a positive factor times) the derivative.]
    // * if this is a inflection point and edge_sample_Lw is above the inflection point,
    //     inflection_point = edge_sample_w - learning_rate * cubic_.derivative(edge_sample_w)
    //     [if we are left/right and above an inflection point, the derivative will be negative/positive;
    //     therefore we need to subtract (a positive factor times) the derivative.]
    // * if this is a inflection point and edge_sample_Lw is below the inflection point,
    //     inflection_point = edge_sample_w + learning_rate * cubic_.derivative(edge_sample_w)
    //     [if we are left/right and below an inflection point, the derivative will be positive/negative;
    //     therefore we need to add (a positive factor times) the derivative.]
    //
    // The point of all this sign toggling is to make sure that learning_rate itself is positive.
    // We can therefore also calculate a `beta` such that:
    //
    //   beta = (edge_sample_w - critical_point_w) / cubic_.derivative(edge_sample_w)
    //
    // and then simply set learning_rate = abs(beta).
    //
    // However, in the case of a maximum/minimum there is a chance that we'd try to divide by zero
    // if edge_sample_w resides in the other extreme (it won't/shouldn't be in the critical_point).
    // Therefore it seems to make more sense to assume a constant for beta, equal to the limit where
    // edge_sample_w approaches critical_point_w: basically assume a parabolic fit in the critical_point_w.
    //
    // Let f(w) = a + b w + c * w^2 be the approximation of the cubic near critical_point_w,
    // where (therefore) critical_point_w = -b / (2c), or b = -2c critical_point_w.
    // Let delta = edge_sample_w - critical_point_w, aka edge_sample_w = critical_point_w + delta.
    // Then beta = delta / (b + 2c (critical_point_w + delta)) =
    //           = delta / (-2c critical_point_w + 2c (critical_point_w + delta)) = 1 / (2c).
    // Where 2c is basically the second derivative in the critical point:
    //
    //   beta = 1 / cubic_.second_derivative(critical_point_w).

    double beta;

    if (type_ == CriticalPointType::inflection_point)
      beta = (members_[side].edge_sample_w - critical_point_w_) / cubic_.derivative(members_[side].edge_sample_w);
    else
      beta = 1 / cubic_.second_derivative(critical_point_w_);

    // Return one tenth of this value.
    return 0.1 * std::abs(beta);
  }

  double or_zero(int side) const { /* Never return really zero */ return std::max(epsilon, value(side)); }

  // Return true if step is significantly smaller than the scale.
  // This is used to determine whether to add a new sample to the history or to replace an existing entry.
  bool negligible(double step) const
  {
    return std::abs(step) < std::max(epsilon, 0.001 * (value(left) + value(right)));
  }

  // Return a value with the same sign as step that is just large enough for negligible to return false;
  double make_significant(double step) const
  {
    double abs_step = 0.00101 * (value(left) + value(right));
    return step < 0.0 ? -abs_step : abs_step;
  }

  // Return true if step is basically zero.
  // This is used to determine if the single sample that we have sits in an extreme,
  // so that a normal gradient descent won't do anything.
  static bool almost_zero(double w, double step)
  {
    return std::abs(step) < epsilon || std::abs(step) < 1e-6 * std::abs(w);
  }

  ScaleUpdate update(std::array<Sample const*, 2> const& relevant_samples, int current_index,
      math::CubicPolynomial const& new_approximation, bool saw_more_than_two_relevant_samples)
  {
    DoutEntering(dc::notice, "Scale::update(" << relevant_samples << ", " << current_index << ", " << new_approximation <<
        ", " << std::boolalpha << saw_more_than_two_relevant_samples << ")");

    // Is this not true?
    ASSERT(saw_more_than_two_relevant_samples || (!members_[left].valid && !members_[right].valid));

    int const prev_sample_index = 1 - current_index;

#if 0
    // Get the x coordinate (w value) of the vertex of the new parabola.
    double const new_v_x = new_parabola.vertex_x();
    Dout(dc::notice, "new_v_x = " << new_v_x);

    // Pick the sample that is horizontally the furthest away from the new vertex.
    Sample const* edge_sample =
      std::abs(relevant_samples[prev_sample_index]->w() - new_v_x) > std::abs(relevant_samples[current_index]->w() - new_v_x) ?
        relevant_samples[prev_sample_index] : relevant_samples[current_index];

    Dout(dc::notice, "The following relevant samples are passed:");
    for (int i = 0; i < 2; ++i)
    {
      Dout(dc::notice|continued_cf, "  " << i << " : " << *relevant_samples[i]);
      if (edge_sample == relevant_samples[i])
        Dout(dc::continued, " (edge sample)");
      Dout(dc::finish, "");
    }

    // If scale is yet unknown, then simply initialize it.
    if (!has_sample_)
    {
      Dout(dc::notice, "Initializing (has_sample_ == false)");

      scale_ = edge_sample->w() - new_v_x;
      Dout(dc::notice, "scale was set to " << *edge_sample << " - " << new_v_x << " = " << scale_);
      edge_sample_w_ = edge_sample->w();
      edge_sample_Lw_ = edge_sample->Lw();
      parabola_ = new_parabola;
      has_sample_ = true;

      return ScaleUpdate::initialized;
    }

    Dout(dc::notice, "The current scale = " << scale_);

    // Get the x coordinate (w value) of the vertex of the old parabola.
    double const old_v_x = parabola_.vertex_x();
    Dout(dc::notice, "old_v_x = " << old_v_x);

    // If the new sample is closer to the (old) vertex then we just jumped to it.
    // Update the scale accordingly.
    if (std::abs(relevant_samples[current_index]->w() - old_v_x) < std::abs(relevant_samples[prev_sample_index]->w() - old_v_x))
    {
      Dout(dc::notice, "The new w (" << *relevant_samples[current_index] << ") is closest to the old vertex (at " << old_v_x << ")");

      // Get the y-coordinate at the w value of the stored edge sample, according to the new parabola.
      double const new_Lw_stored_edge_sample = new_parabola(edge_sample_w_);
      // Get the vertical distance from the stored edge sample to the new parabola.
      double const distance_stored_edge_sample_to_parabola = std::abs(new_Lw_stored_edge_sample - edge_sample_Lw_);
      // Get the vertical distance from the stored edge sample to the vertex of the new parabola
      // (equal to abs(new_parabola.height(edge_sample_w_))).
      double const abs_stored_edge_sample_height = std::abs(new_Lw_stored_edge_sample - new_parabola.vertex_y());
      // If the stored sample vertically deviates more than 10%, we need to adjust the stored sample values.
      if (distance_stored_edge_sample_to_parabola > 0.1 * abs_stored_edge_sample_height)
      {
        Dout(dc::notice, "Discarded old sample at L(" << edge_sample_w_ << ") = " << edge_sample_Lw_ << " because " <<
            distance_stored_edge_sample_to_parabola << " > 0.1 * " << abs_stored_edge_sample_height);

        static constexpr int left = -1;
        static constexpr int right = 1;
        // On which side of the vertex is the old sample?
        int hside = edge_sample_w_ < new_v_x ? left : right;

        int count;
        std::array<double, 4> toggles;
        bool close = new_parabola.equal_intervals(parabola_, toggles, count);
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
          new_edge_w = std::min(close ? edge_sample_w_ : toggles[i], edge_sample->w());
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
          new_edge_w = std::max(close ? edge_sample_w_ : toggles[i], edge_sample->w());
        }

        edge_sample_w_ = new_edge_w;
        edge_sample_Lw_ = utils::almost_equal(new_edge_w, edge_sample->w(), 1e-3) ? edge_sample->Lw() : new_parabola(new_edge_w);
      }
      scale_ = edge_sample_w_ - new_v_x;
      Dout(dc::notice, "scale was set to " << edge_sample_w_ << " - " << new_v_x << " = " << scale_);

      parabola_ = new_parabola;

      return ScaleUpdate::towards_vertex;
    }

    Sample const* current = relevant_samples[current_index];
    // Calculate the potential new scale relative to the vertex of the canonical parabola.
    double new_scale = current->w() - old_v_x;
    // If the scale could grow, see if the new sample is close enough to the parabola.
    if (std::abs(new_scale) > std::abs(scale_))
    {
      // Get the y-coordinate at the w value of the new sample, according to the parabola.
      double const parabola_Lw_at_current = parabola_(current->w());
      // Get the vertical distance from the new sample to the new parabola.
      double const distance_current_to_parabola = std::abs(parabola_Lw_at_current - current->Lw());
      // Get the vertical distance from the new sample to the vertex of the parabola.
      double const abs_current_height = std::abs(current->Lw() - parabola_.vertex_y());
      // If the new sample vertically deviates more than 10%, we can't use the new sample as new edge sample.
      if (distance_current_to_parabola > 0.1 * abs_current_height)
        return ScaleUpdate::disconnected;
      // Get the derivative at the w value of the new sample, according to the parabola and adjust the scale accordingly.
      double const parabola_dLdw_at_current = parabola_.derivative(current->w());
      // Get the real derivative at this w value.
      double const dLdw_at_current = current->dLdw();
      // The sin of the angle between the tangent lines at this point (with slopes corresponding to these derivatives)
      // is given by: sin(delta_angle) = (s_1 - s_2) / sqrt((1 + (s_1)^2)(1 + (s_2)^2)).
      // Calculate the square of the sinus of the angle between these lines.
      double const s1s2 = parabola_dLdw_at_current * dLdw_at_current;
      double const s12 = parabola_dLdw_at_current * parabola_dLdw_at_current;
      double const s22 = dLdw_at_current * dLdw_at_current;
      double const sin2 = (s12 - 2.0 * s1s2 + s22) / ((1.0 + s12) * (1.0 + s22));
      if (sin2 < 0.03)
      {
        scale_ = new_scale;
        edge_sample_w_ = current->w();
        edge_sample_Lw_ = current->Lw();
        Dout(dc::notice, "scale was set to " << *current << " - " << old_v_x << " = " << scale_);
      }
    }
#else
    // Implement
    ASSERT(false);
#endif
    return ScaleUpdate::away_from_vertex;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{scale:[" << members_[left].scale << ", " << members_[right].scale << "], has_sample:[" <<
      std::boolalpha << members_[left].valid << ", " << members_[right].valid << "]}";
  }
#endif
};

} // namespace gradient_descent
