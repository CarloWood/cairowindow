#pragma once

#include "Sample.h"
#include "HorizontalDirection.h"
#include "ExtremeType.h"
#include "../CubicPolynomial.h"
#include "utils/debug_ostream_operators.h"
#include "utils/almost_equal.h"
#include "utils/macros.h"
#include "utils/DEVector.h"
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
  towards_cp,                   // When we already had two samples (and therefore an "old" cubic approximation)
                                // and the new sample is at the critical point of the new cubic.
  away_from_cp,                 // When we already had two samples (and therefore an "old" cubic approximation)
                                // and the new sample was moved away from the critical point.
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
  utils::DEVector<SampleData> samples_; // All samples previously passed to add.

 public:
  Scale() = default;
  Scale(Scale const&) = delete;
  Scale(Scale&& orig) : valid_{orig.valid_}, value_{orig.value_}, left_edge_{orig.left_edge_}, right_edge_{orig.right_edge_},
    cubic_{std::move(orig.cubic_)}, type_{orig.type_}, critical_point_w_{orig.critical_point_w_}, samples_{std::move(orig.samples_)} { }

 public:
  Scale& operator=(Scale const& scale) = delete;

  // Accessors.
  bool is_valid() const { return valid_; }
  double value() const { ASSERT(valid_); return value_; }
  double or_zero() const { return valid_ ? value_ : epsilon; }
  math::CubicPolynomial const& cubic() const { return cubic_; }
  CriticalPointType type() const { return type_; }
  double inflection_point_w() const { return cubic_.inflection_point(); }
  double critical_point_w() const { return critical_point_w_; }
  double left_edge_sample_w() const { ASSERT(valid_); return samples_[left_edge_].w; }
  double right_edge_sample_w() const { ASSERT(valid_); return samples_[right_edge_].w; }
  double left_edge_sample_Lw() const { ASSERT(valid_); return samples_[left_edge_].sample->Lw(); }
  double right_edge_sample_Lw() const { ASSERT(valid_); return samples_[right_edge_].sample->Lw(); }

  Sample const* get_nearest_sample(Sample const* current) const
  {
    // This function is called on a Scale that is part of a LocalExtreme. It can not be empty.
    ASSERT(samples_.size() > 2);        // Two samples to make a cubic, plus the disconnected current sample.
    auto iter = std::find_if(samples_.begin(), samples_.end(), [current](SampleData const& data){ return data.sample == current; });
    // current was already added to this Scale.
    ASSERT(iter != samples_.end());

    utils::DEVector<SampleData>::const_iterator right_iter = iter + 1;
    if (iter == samples_.begin())
      return right_iter->sample;
    utils::DEVector<SampleData>::const_iterator left_iter = iter - 1;
    if (right_iter == samples_.end())
      return left_iter->sample;

    return ((iter->w - left_iter->w < right_iter->w - iter->w) ? left_iter : right_iter)->sample;
  }

  // Return a directional scale: minus the scale if direction is `left` and plus the scale if direction is `right`.
  double step(HorizontalDirectionToInt direction) const
  {
    DoutEntering(dc::notice, "Scale::step(" << direction << ")");
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
      Dout(dc::notice, "Setting direction to " << dir << " as the furthest sample from " << critical_point_w_ <<
          " is " << *samples_[far_edge].sample);
    }
    return dir * value_;
  }

  void reset()
  {
    DoutEntering(dc::notice, "Scale::reset()");
    samples_.clear();
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

 public:
  static constexpr int Li = 0;
  static constexpr int Ci = 1;
  static constexpr int Ri = 2;
  static constexpr int Ii = 3;
  using LCRI_type = std::array<double, 4>;      // An array with the L, C, R and I values.

  static constexpr unsigned char const classify_map[13] = {4, 4, 1, 4, 3, 2, 0, 3, 4, 1, 4, 4};
  static int classify(LCRI_type const& LCRI)
  {
    return classify_map[(LCRI[Ii] < LCRI[Li]) + (LCRI[Ii] < LCRI[Ri]) +
      4 * ((LCRI[Ci] < LCRI[Li]) + (LCRI[Ci] < LCRI[Ri])) + (LCRI[Ii] < LCRI[Ci])];
  }

 private:
  // Return the weighted average distance of L and R to C.
  static double weighted_average(LCRI_type const& LCRI)
  {
    ASSERT(LCRI[Li] < LCRI[Ri]);        // In fact, LCRI[Li] < LCRI[Ci] < LCRI[Ri].
    DoutEntering(dc::notice, "Scale::weighted_average(" << LCRI[Li] << ", " << LCRI[Ci] << ", " << LCRI[Ri] << ")");
    return (utils::square(LCRI[Ci] - LCRI[Li]) + utils::square(LCRI[Ri] - LCRI[Ci])) / (LCRI[Ri] - LCRI[Li]);
  }

  double calculate_value() const
  {
    // The (value of the) scale is losely defined at the distance over which w can change
    // without L(w) starting to deviate significantly from the (current) approximation / cubic.
    //
    // Moreover, because of the nature of the algorithm, which approaches a critical point
    // from one side, we don't want the scale to be much different if one side hasn't be
    // explored yet.
    //
    // For example,
    //
    //      cp          : critical point
    //     l|           : left-most sample
    //     ↓↓
    //     .-. r        : right-most sample
    //    /   \↓
    //   /     \
    //
    // If l is very close to cp, or non-existant, which gives a distance of 0, then
    // we just want scale to be `r - cp`; assuming that there is no reason that
    // L(w) would look different on the left side of cp if we never explored that side.
    //
    // But if `cp - l` is significant, although less than `r - cp` the chance increases
    // that we DID explore the left side and there is too much deviation left of l.
    //
    // Therefore, we return a "weighted" average between `cp - l` and `r - cp`, weighted
    // with their own value: (w1*A + w2*B) / (w1 + w2), where we'll use w1=A and w2=B.
    // In this case A = cp - l and B = r - cp. Hence,
    //
    //   value = (A^2 + B^2) / (A + B)
    //
    // Where A and B are positive if l and r are on opposite sides of cp.
    //
    // However, if l and r are on the same side of cp, for example,
    //
    //      cp          : critical point
    //     r|           : right-most sample (now B is negative)
    //     ↓↓
    //   l .-.          : left-most sample
    //   ↓/   \
    //   /     \
    //
    // Then we want the returned scale to be `r - l`.
    //
    // Extension with inflection point
    // -------------------------------
    //
    // Let C be the critical point at critical_point_w_.
    // Let I be the inflection point of the cubic.
    // Then the following situations are possible:
    //
    //                 C/I                                    : C and I on top of eachother (C is the inflection point)
    //                                                          this happens when the cubic has no local extremes.
    //
    //           I     C                                      : The inflection point is on the left of the critical point.
    //                                                          \      C                    /
    //                                                           \    /\             /\    /
    //                                                            \  I  \     or    /  I  /
    //                                                             \/    \         /    \/
    //                                                                    \       /      C
    //                 C     I                                : The inflection point is on the right of the critical point.
    //                                                          \                    C      /
    //                                                           \    /\             /\    /
    //                                                            \  I  \     or    /  I  /
    //                                                             \/    \         /    \/
    //                                                             C      \       /
    //
    // Let the left_edge_ sample be l and the right_edge_ sample be r.
    // Placing l and r in each of the above cases we have the possibilities:
    //
    //                                      Class (value returned by classify(l, I, C, r)).
    //   l < r < C = I   : r - l                4
    //   C = I < l < r   : r - l                4
    //   l < r < I < C   : r - l                4
    //   I < l < r < C   : r - l                4
    //   I < C < l < r   : r - l                4
    //   l < r < C < I   : r - l                4
    //   C < l < r < I   : r - l                4
    //   C < I < l < r   : r - l                4
    //   l < I < r < C   : wa(l, I, r)          1
    //   C < l < I < r   : wa(l, I, r)          1
    //   l < C = I < r   : wa(l, C, r)          3 - however, classify returns 2.
    //   I < l < C < r   : wa(l, C, r)          3
    //   l < C < r < I   : wa(l, C, r)          3
    //   l < I < C < r   : wa(I, C, r)          0
    //   l < C < I < r   : wa(l, C, I)          2
    //
    // where wa(l, C, r) is the weighted average around C: ((C - l)^2 + (r - C)^2) / (r - l) etc.
    //
    // Note that the Class values have been chosen such that we can simply replace element `Class`
    // in the LCRI array with the inflection point, replacing respectively l, C or r (or I with itself
    // in the case of Class 3) and then call weighted_average with that altered array.

    DoutEntering(dc::notice, "Scale::calculate_value()");
    double const inflection_point_w = cubic_.inflection_point();
    Dout(dc::notice, "left_edge_w = " << samples_[left_edge_].w <<
        "; critical_point_w = " << critical_point_w_ <<
        "; right_edge_w = " << samples_[right_edge_].w <<
        "; inflection_point_w = " << inflection_point_w);
    // critical_point_w_ should be initialized.
    ASSERT(type_ != CriticalPointType::none);
    LCRI_type LCRI = {samples_[left_edge_].w, critical_point_w_, samples_[right_edge_].w, inflection_point_w};
    int Class = classify(LCRI);
    double result;
    if (Class == 4)
      result = LCRI[Ri] - LCRI[Li];
    else
    {
      if (AI_LIKELY(type_ != CriticalPointType::inflection_point))
        LCRI[Class] = LCRI[Ii];
      result = weighted_average(LCRI);
    }
    Dout(dc::notice, "Returning value: " << result);
    return result;
  }

 public:
  ScaleUpdate update(ExtremeType extreme_type, std::array<Sample const*, 2> const& relevant_samples, int current_index,
      math::CubicPolynomial const& new_cubic, bool saw_more_than_two_relevant_samples, bool local_extreme)
  {
    // Note: if this an update of an already found local extreme, then `extreme_type` is
    // the current extreme type (that matches the already found local extreme).
    // Otherwise it is the next_extreme_type: the extreme type that we're looking for.
#ifdef CWDEBUG
    DoutEntering(dc::notice|continued_cf, "Scale::update(" << extreme_type << ", {");
    char const* separator = "";
    for (int i = 0; i < relevant_samples.size(); ++i)
    {
      Dout(dc::continued, separator << '&' << *relevant_samples[i]);
      if (i == current_index)
        Dout(dc::continued, " (current)");
      separator = ", ";
    }
    Dout(dc::finish, "}, " << current_index << ", " << new_cubic << ", " << std::boolalpha << saw_more_than_two_relevant_samples <<
        ", " << local_extreme << ")");
#endif

    // If we didn't have more than two relevant samples, then these two are the first two and this scale should still be invalid.
    ASSERT(saw_more_than_two_relevant_samples || !valid_);

    int const prev_sample_index = 1 - current_index;

    if (type_ == CriticalPointType::none)
    {
      // In the case of a local_extreme, type_ is type of that extreme, so it can't be 'none'.
      ASSERT(!local_extreme);
      // extreme_type was set by find_extreme.
      if (extreme_type == ExtremeType::unknown)
      {
        //FIXME: is this correct?
        type_ = CriticalPointType::inflection_point;
      }
      else
        type_ = extreme_type == ExtremeType::minimum ? CriticalPointType::minimum : CriticalPointType::maximum;
    }
    else
    {
      // The extreme types of this Scale and whatever Algorithm is looking for should correspond.
      // Can this fail because one is unknown, or inflection_point? If so, what to do here?
      ASSERT((extreme_type == ExtremeType::minimum) == (type_ == CriticalPointType::minimum));
      ASSERT((extreme_type == ExtremeType::maximum) == (type_ == CriticalPointType::maximum));
    }

    [[maybe_unused]] double new_cp_w;
    if (!local_extreme)
    {
      // Get the w coordinate of the critical point of the new cubic.
      std::array<double, 2> extremes;
      int number_of_extremes = new_cubic.get_extremes(extremes, false);
      cubic_ = new_cubic;
      new_cp_w = (number_of_extremes == 2 && type_ == CriticalPointType::maximum) ? extremes[1] : extremes[0];
      Dout(dc::notice, "new_cp_w = " << new_cp_w);
    }
    else
    {
      // The approximation wasn't changed if this is a local extreme.
      ASSERT(cubic_ == new_cubic);
    }

    if (!saw_more_than_two_relevant_samples)
    {
      // If we already found a local extreme then we had already two samples before this one.
      ASSERT(!local_extreme);
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

    Sample const* current = relevant_samples[current_index];
    double w = current->w();
#ifdef CWDEBUG
    {
      double prev_w = relevant_samples[1 - current_index]->w();
      bool towards_cp;
      if (local_extreme)
        towards_cp = std::abs(w - critical_point_w_) < std::abs(prev_w - critical_point_w_);
      else
        towards_cp = std::abs(w - new_cp_w) < std::abs(prev_w - new_cp_w);
      // Shouldn't we always move towards the critical point when this isn't a found local extreme,
      // and away from it when it is already found?
      //FIXME: uncomment this.
      //ASSERT(towards_cp != local_extreme);
    }
#endif

    bool at_back = samples_.back().w <= w;
    bool at_front = samples_.front().w > w;
    std::size_t index;          // The index into samples_ of current, after current was inserted.
    if (at_front)
    {
      index = 0;
      samples_.emplace_front(w, current);
      ++left_edge_;
      ++right_edge_;
    }
    else if (at_back)
    {
      index = samples_.size();
      samples_.emplace_back(w, current);
    }
    else
    {
      // Support w = 2, and we have samples [-1, 1, 3, 4]
      // then iter will point to                    ^--   the first element larger than w.
      // We then insert before that element to get:
      // [-1, 1, 2, 3, 4] - keeping the vector sorted.
      //   0  1  2 <-- index.
      auto iter = std::find_if(samples_.begin(), samples_.end(), [w](SampleData const& element) { return element.w > w; });
      index = std::distance(samples_.begin(), iter);
      samples_.emplace(iter, w, current);
#if 0
      if (index <= left_edge_)
      {
        ++left_edge_;
        ++right_edge_;
      }
      else if (index <= right_edge_)
        ++right_edge_;
#else
      // We shouldn't try a sample in the same interval that was already tested in this case.
      // Aka, if this is a local extreme then at_back or at_front is expected to be true.
//FIXME: uncomment and fix?
//      ASSERT(!local_extreme);
#endif
    }

    double cp_Lw = cubic_(critical_point_w_);
    if (local_extreme)
    {
      double Lw = current->Lw();
      double Aw = cubic_(w);
      if (std::abs(Aw - Lw) > 0.1 * std::abs(Lw - cp_Lw))
      {
        Dout(dc::notice, "Ignoring new sample " << *current << " because cubic_(" << w << ") = " <<
            Aw << " and |" << Aw << " - " << Lw << "| = " << (std::abs(Aw - Lw)) << " > 0.1 * |" << Lw << " - " << cp_Lw << "| = " <<
            0.1 * std::abs(Lw - cp_Lw));
        return ScaleUpdate::disconnected;
      }
      else
      {
        // It seems expected that current is right next to another matching, already added, sample, no?
//FIXME: uncomment and fix?
//        ASSERT(index == left_edge_ - 1 || index == right_edge_ + 1);
        if (index < left_edge_)
          --left_edge_;
        else
          ++right_edge_;
      }
    }
    else
    {
      // Find the first element in samples_ that has a w value that is larger than current and then insert current in front of that.
      critical_point_w_ = new_cp_w;
      // The current approximation cubic is based on current.
      ASSERT(utils::almost_equal(cubic_(w), current->Lw(), 10e-5));
      left_edge_ = right_edge_ = index;
      while (left_edge_ > 0)
      {
        --left_edge_;
        double Lw = samples_[left_edge_].sample->Lw();
        double Aw = cubic_(samples_[left_edge_].w);
        if (std::abs(Aw - Lw) > 0.1 * std::abs(Lw - cp_Lw))
        {
          Dout(dc::notice, "Ignoring sample " << *samples_[left_edge_].sample << " because cubic_(" << samples_[left_edge_].w << ") = " <<
              Aw << " and |" << Aw << " - " << Lw << "| = " << (std::abs(Aw - Lw)) << " > 0.1 * |" << Lw << " - " << cp_Lw << "| = " <<
              0.1 * std::abs(Lw - cp_Lw));
          ++left_edge_;
          break;
        }
        Dout(dc::notice, "Including sample " << *samples_[left_edge_].sample << " because cubic_(" << samples_[left_edge_].w << ") = " <<
            Aw << " and |" << Aw << " - " << Lw << "| = " << (std::abs(Aw - Lw)) << " <= 0.1 * |" << Lw << " - " << cp_Lw << "| = " <<
            0.1 * std::abs(Lw - cp_Lw));
      }
      Dout(dc::notice, "left_edge_ is now " << left_edge_);
      while (right_edge_ + 1 < samples_.size())
      {
        ++right_edge_;
        double Lw = samples_[right_edge_].sample->Lw();
        double Aw = cubic_(samples_[right_edge_].w);
        if (std::abs(Aw - Lw) > 0.1 * std::abs(Lw - cp_Lw))
        {
          Dout(dc::notice, "Ignoring sample " << *samples_[right_edge_].sample << " because cubic_(" << samples_[right_edge_].w << ") = " <<
              Aw << " and |" << Aw << " - " << Lw << "| = " << (std::abs(Aw - Lw)) << " > 0.1 * |" << Lw << " - " << cp_Lw << "| = " <<
              0.1 * std::abs(Lw - cp_Lw));
          --right_edge_;
          break;
        }
        Dout(dc::notice, "Including sample " << *samples_[right_edge_].sample << " because cubic_(" << samples_[right_edge_].w << ") = " <<
            Aw << " and |" << Aw << " - " << Lw << "| = " << (std::abs(Aw - Lw)) << " <= 0.1 * |" << Lw << " - " << cp_Lw << "| = " <<
            0.1 * std::abs(Lw - cp_Lw));
      }
      Dout(dc::notice, "right_edge_ is now " << right_edge_);
    }

    value_ = calculate_value();
    return local_extreme ? ScaleUpdate::away_from_cp : ScaleUpdate::towards_cp;
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
