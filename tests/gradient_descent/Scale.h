#pragma once

#include "Sample.h"
#include "HorizontalDirection.h"
#include "ExtremeType.h"
#include "../CubicPolynomial.h"
#include "utils/debug_ostream_operators.h"
#include "utils/almost_equal.h"
#include "utils/macros.h"
#include <boost/container/devector.hpp>
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
  towards_cp,                   // When we already had two samples (and therefore an "old" parabolic approximation)
                                // and the new sample is at the vertex of the new parabola.
  away_from_cp,                 // When we already had two samples (and therefore an "old" parabolic approximation)
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
    AI_CASE_RETURN(ScaleUpdate::towards_cp);
    AI_CASE_RETURN(ScaleUpdate::away_from_cp);
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
  none,
  minimum,
  maximum,
  inflection_point
};

#ifdef CWDEBUG
inline std::string to_string(CriticalPointType critical_point_type)
{
  switch (critical_point_type)
  {
    AI_CASE_RETURN(CriticalPointType::none);
    AI_CASE_RETURN(CriticalPointType::minimum);
    AI_CASE_RETURN(CriticalPointType::maximum);
    AI_CASE_RETURN(CriticalPointType::inflection_point);
  }
  AI_NEVER_REACHED
}

inline std::ostream& operator<<(std::ostream& os, CriticalPointType critical_point_type)
{
  return os << to_string(critical_point_type);
}
#endif

class Scale
{
 public:
  static constexpr double epsilon = 1e-30;

  struct SampleData {
    double w;
    Sample const* sample;
  };

 protected:
  bool valid_{false};                   // True iff the members below are valid.
  double value_{};                      // Cache of value(). An indication of what changes to w are significant,
                                        // relative to the critical point. This value is always larger than zero.
  int left_edge_;                       // Index into samples_ of the left-most sample that still matches the cubic approximation.
  int right_edge_;                      // Index into samples_ of the right-most sample that still matches the cubic approximation.

  math::CubicPolynomial cubic_;         // The last (previous) third degree polynomial fit (passed to update).
  CriticalPointType type_{CriticalPointType::none};     // Whether distances are measured relative to the minimum, maximum or
                                        // the inflection point (if the cubic has no extremes).
  double critical_point_w_;             // Cached value of the x-coordinate of the critical point.
  boost::container::devector<SampleData> samples_;      // All samples previously passed to add.

 public:
  Scale() = default;

 public:
  Scale& operator=(Scale const& scale) = delete;

  // Accessors.
  bool is_valid() const { return valid_; }
  double value() const { ASSERT(valid_); return value_; }
  double or_zero() const { return valid_ ? value_ : epsilon; }
  math::CubicPolynomial const& cubic() const { return cubic_; }
  CriticalPointType type() const { return type_; }
  double critical_point_w() const { return critical_point_w_; }
  double left_edge_sample_w() const { ASSERT(valid_); return samples_[left_edge_].w; }
  double right_edge_sample_w() const { ASSERT(valid_); return samples_[right_edge_].w; }
  double left_edge_sample_Lw() const { ASSERT(valid_); return samples_[left_edge_].sample->Lw(); }
  double right_edge_sample_Lw() const { ASSERT(valid_); return samples_[right_edge_].sample->Lw(); }

  // Return a directional scale: minus the scale if direction is `left` and plus the scale if direction is `right`.
  double step(HorizontalDirectionToInt direction) const
  {
    ASSERT(valid_);
    int dir = direction;
    if (dir == 0)
    {
      // If the direction passed is undecided, then set dir to the direction from the sample that is
      // the furthest away from the critical_point towards the critical point (even including samples
      // that don't fit the current cubic).
      int far_edge = std::abs(samples_.front().w - critical_point_w_) > std::abs(samples_.back().w - critical_point_w_) ?
        0 : samples_.size() - 1;
      dir = std::copysign(1.0, critical_point_w_ - samples_[far_edge].w);
    }
    return dir * value_;
  }

  void reset()
  {
    DoutEntering(dc::notice, "Scale::reset()");
    type_ = CriticalPointType::none;
    valid_ = false;
  }

  // Return true if step is significantly smaller than the scale.
  // This is used to determine whether to add a new sample to the history or to replace an existing entry.
  bool negligible(double step) const
  {
    return std::abs(step) < std::max(epsilon, 0.001 * or_zero());
  }

  // Return a value with the same sign as step that is just large enough for negligible to return false;
  double make_significant(double step) const
  {
    double abs_step = 0.00101 * value();
    return step < 0.0 ? -abs_step : abs_step;
  }

  // Return true if step is basically zero.
  // This is used to determine if the single sample that we have sits in an extreme,
  // so that a normal gradient descent won't do anything.
  static bool almost_zero(double w, double step)
  {
    return std::abs(step) < epsilon || std::abs(step) < 1e-6 * std::abs(w);
  }

  double calculate_value() const
  {
    DoutEntering(dc::notice, "Scale::calculate_value()");
    Dout(dc::notice, "left_edge_w = " << samples_[left_edge_].w << "; right_edge_w = " << samples_[right_edge_].w <<
        "; critical_point_w = " << critical_point_w_);
    // critical_point_w_ should be initialized.
    ASSERT(type_ != CriticalPointType::none);
    double L = samples_[left_edge_].w - critical_point_w_;
    double R = samples_[right_edge_].w - critical_point_w_;
    Dout(dc::notice, "L = " << L << "; R = " << R);
    double sign_L = std::copysign(1.0, L);
    double sign_R = std::copysign(1.0, R);
    double result;
    if (sign_L == sign_R)
      result = sign_L == -1.0 ? -L : R;
    else
    {
      // L must be left of R.
      ASSERT(sign_L == -1.0 && sign_R == 1.0);
      result = (L * L + R * R) / (-L + R);
    }
    Dout(dc::notice, "Returning value: " << result);
    return result;
  }

  ScaleUpdate update(ExtremeType next_extreme_type, std::array<Sample const*, 2> const& relevant_samples, int current_index,
      math::CubicPolynomial const& new_cubic, bool saw_more_than_two_relevant_samples)
  {
#ifdef CWDEBUG
    DoutEntering(dc::notice|continued_cf, "Scale::update(" << next_extreme_type << ", {");
    char const* separator = "";
    for (int i = 0; i < relevant_samples.size(); ++i)
    {
      Dout(dc::continued, separator << '&' << *relevant_samples[i]);
      if (i == current_index)
        Dout(dc::continued, " (current)");
      separator = ", ";
    }
    Dout(dc::finish, "}, " << current_index << ", " << new_cubic << ", " << std::boolalpha << saw_more_than_two_relevant_samples << ")");
#endif

    // If we didn't have more than two relevant samples, then these two are the first two and this scale should still be invalid.
    ASSERT(saw_more_than_two_relevant_samples || !valid_);

    int const prev_sample_index = 1 - current_index;

    if (type_ == CriticalPointType::none)
    {
      // next_extreme_type was set by find_extreme.
      if (next_extreme_type == ExtremeType::unknown)
        //FIXME: is this correct?
        type_ = CriticalPointType::inflection_point;
      else
        type_ = next_extreme_type == ExtremeType::minimum ? CriticalPointType::minimum : CriticalPointType::maximum;
    }
    else
    {
      // The extreme types of this Scale and whatever Algorithm is looking for should correspond.
      // Can this fail because one is unknown, or inflection_point? If so, what to do here?
      ASSERT((next_extreme_type == ExtremeType::minimum) == (type_ == CriticalPointType::minimum));
      ASSERT((next_extreme_type == ExtremeType::maximum) == (type_ == CriticalPointType::maximum));
    }

    // Get the w coordinate of the critical point of the new cubic.
    std::array<double, 2> extremes;
    int number_of_extremes = new_cubic.get_extremes(extremes, false);
    cubic_ = new_cubic;
    double new_cp_w = (number_of_extremes == 2 && type_ == CriticalPointType::maximum) ? extremes[1] : extremes[0];
    Dout(dc::notice, "new_cp_w = " << new_cp_w);

    if (!saw_more_than_two_relevant_samples)
    {
      critical_point_w_ = new_cp_w;
      // This should be the first time this function is called.
      ASSERT(samples_.empty());
      // Add the two samples in order of w value.
      int left_index = relevant_samples[0]->w() < relevant_samples[1]->w() ? 0 : 1;
      int right_index = 1 - left_index;
      samples_.emplace_front(relevant_samples[left_index]->w(), relevant_samples[left_index]);
      samples_.emplace_back(relevant_samples[right_index]->w(), relevant_samples[right_index]);
      left_edge_ = 0;
      right_edge_ = 1;
      value_ = calculate_value();
      valid_ = true;
      return ScaleUpdate::initialized;
    }

    // Find the first element in samples_ that has a w value that is larger than current and then
    // insert current in front of that.
    Sample const* current = relevant_samples[current_index];
    double w = current->w();
    double prev_w = relevant_samples[1 - current_index]->w();
    bool towards_cp = std::abs(w - critical_point_w_) < std::abs(prev_w - critical_point_w_);
    critical_point_w_ = new_cp_w;
    bool at_back = samples_.back().w <= w;
    bool at_front = samples_.front().w > w;
    ++right_edge_;
    if (at_front)
      samples_.emplace_front(w, current);
    else if (at_back)
      samples_.emplace_back(w, current);
    else
    {
      auto iter = std::find_if(samples_.begin(), samples_.end(), [w](SampleData const& element) { return element.w > w; });
      samples_.emplace(iter, w, current);
    }
    value_ = calculate_value();

    return towards_cp ? ScaleUpdate::towards_cp : ScaleUpdate::away_from_cp;
#if 0
    //...

    // Get the w coordinate of the critical point of the old cubic.
    number_of_extremes = cubic_.get_extremes(extremes, false);
    double const old_cp_w = (number_of_extremes == 2 && type_ == CriticalPointType::maximum) ? extremes[1] : extremes[0];
    Dout(dc::notice, "old_cp_w = " << old_cp_w);

    // If the new sample is closer to the (old) vertex then we just jumped to it.
    // Update the scale accordingly.
    if (std::abs(relevant_samples[current_index]->w() - old_cp_w) < std::abs(relevant_samples[prev_sample_index]->w() - old_cp_w))
    {
      Dout(dc::notice, "The new w (" << *relevant_samples[current_index] << ") is closest to the old vertex (at " << old_cp_w << ")");
#endif
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{valid_:" << std::boolalpha << valid_ << ", samples_:{";
    char const* separator = "";
    for (auto&& sample_data : samples_)
    {
      os << separator << *sample_data.sample;
      separator = ", ";
    }
    os << "}";
    if (valid_)
      os << ", left_edge:" << left_edge_ << ", right_edge:" << right_edge_ << ", cubic:" << cubic_;
    if (type_ != CriticalPointType::none)
      os << ", critical_point_w:" << critical_point_w_;
    os << "}";
  }
#endif
};

} // namespace gradient_descent
