#pragma once

#include "Sample.h"
#include "Scale.h"
#include "CubicToNextSampleType.h"
#include "HorizontalDirection.h"
#include "ExtremeType.h"
#include "CriticalPointType.h"
#include "../CubicPolynomial.h"
#include <memory>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include "events/Events.h"
#include <iostream>
#endif

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
struct AlgorithmEventType;
#endif

// SampleNode's form a linked list, stored in Algorithm, of all Sample's done
// so far. They are strictly ordered by Sample::w_ value: the w value of prev_
// is less than, and the w value of next_ is larger than that of this sample.
//
// Any additional data stored in SampleNode refers to the cubic that fits
// through this Sample and the next_ Sample.
class SampleNode : public Sample
{
 public:
  static constexpr int explored_left = 1;               // Bit mask used for explored_.
  static constexpr int explored_right = 2;

  using list_type = std::list<SampleNode>;              // The actual list is member of ExtremeChain.
  using iterator = list_type::iterator;
  using const_iterator = list_type::const_iterator;

#if CW_DEBUG
 public:
  // Return value of find_extreme when no extreme is available.
  static constexpr double uninitialized_magic = 12345678.876543211;
#endif

 private:
  mutable math::CubicPolynomial cubic_;                         // The cubic that fits this and the next Sample.
  mutable CubicToNextSampleType type_;                          // The type of this cubic.
  mutable Scale scale_;                                         // The scale that belongs to this cubic.
  mutable const_iterator next_sample_;                          // The right-sample used for cubic_, copy of what was passed to initialize_cubic.
  mutable ExtremeType local_extreme_{ExtremeType::unknown};     // Set to minimum or maximum if this is a local extreme. Otherwise set to unknown.
  mutable double extreme_Lw_;                                   // The Lw coordinate of the local extreme (w is stored in the scale_).
  mutable int explored_{0};                                     // Bit mask 1: exploration to the left of this extreme has started.
                                                                // Bit mask 2: same, on the right.
                                                                // Only valid if this is a local extreme.
  mutable double opposite_direction_w_;                         // The direction to explore if we jump back here as best minimum, only valid if
                                                                // local_extreme_ isn't unknown.

 public:
  SampleNode(Sample&& sample) : Sample(std::move(sample)), type_(CubicToNextSampleType::unknown) { }

  void initialize_cubic(const_iterator next
      COMMA_CWDEBUG_ONLY(events::Server<AlgorithmEventType>& event_server, bool this_is_last)) const;

 public:
  double find_extreme(Sample const& next, ExtremeType& extreme_type) const;

  void set_scale(CriticalPointType type, double critical_point_w, double left_edge_w, double right_edge_w) const
  {
    scale_.set(type, critical_point_w, left_edge_w, right_edge_w, cubic_.inflection_point());
  }

  const_iterator next_sample() const
  {
    // Call initialize_cubic before using this member function.
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return next_sample_;
  }

  //---------------------------------------------------------------------------
  // Accessors for the cubic.

  // Returns the cubic that was fitted between this and the next sample.
  math::CubicPolynomial const& cubic() const { return cubic_; }

  // The scale of this cubic.
  Scale const& scale() const { return scale_; }

  // Returns the type of the cubic that was fitted between this and the next sample.
  CubicToNextSampleType type() const { return type_; }

  // Returns the w value of the underlying extreme.
  double extreme_w() const
  {
    // Call initialize_cubic before using this member function.
    ASSERT(type_ != CubicToNextSampleType::unknown);
    // Call set_local_extreme before calling this function.
    ASSERT(local_extreme_ != ExtremeType::unknown);
    return scale_.critical_point_w();
  }

  // Returns true if this type falls in the 'rising' class.
  bool is_rising() const
  {
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return (static_cast<int>(type_) & rising_bit) != 0 ||
      ((static_cast<int>(type_) & (minimum_bit|maximum_bit)) != 0 && next_sample_->Lw() > Lw());
  }

  // Returns true if this type falls in the 'falling' class.
  bool is_falling() const
  {
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return (static_cast<int>(type_) & falling_bit) != 0 ||
      ((static_cast<int>(type_) & (minimum_bit|maximum_bit)) != 0 && next_sample_->Lw() < Lw());
  }

  // Returns true if the cubic fitted between this and the next sample has an extreme.
  bool has_extreme(ExtremeType extreme_type) const
  {
    ASSERT(extreme_type != ExtremeType::unknown);
    return extreme_type == ExtremeType::minimum ? has_minimum(type_) : has_maximum(type_);
  }

  // Scale estimates that can be used when scale is not available yet.
  double w_scale_estimate() const { ASSERT(type_ != CubicToNextSampleType::unknown); ASSERT(next_sample_->w() > w()); return next_sample_->w() - w(); }
  double Lw_scale_estimate() const { return std::abs(next_sample_->Lw() - Lw()); }
  double dLdw_scale_estimate() const { return std::abs((next_sample_->Lw() - Lw()) / (next_sample_->w() - w())); }

  bool has_left_min() const
  {
    return (static_cast<int>(type_) & left_min_bit) != 0;
  }

  bool has_left_max() const
  {
    return (static_cast<int>(type_) & left_max_bit) != 0;
  }

  bool has_right_min() const
  {
    return (static_cast<int>(type_) & right_min_bit) != 0;
  }

  bool has_right_max() const
  {
    return (static_cast<int>(type_) & right_max_bit) != 0;
  }

  //---------------------------------------------------------------------------
  // Local extreme functions.

  void set_local_extreme(ExtremeType local_extreme) const
  {
    ASSERT(scale_.type() != CriticalPointType::none);
    extreme_Lw_ = cubic_(scale_.critical_point_w());
    local_extreme_ = local_extreme;
  }

  bool is_local_extreme() const { return local_extreme_ != ExtremeType::unknown; }
  ExtremeType get_extreme_type() const { ASSERT(local_extreme_ != ExtremeType::unknown); return local_extreme_; }
  double extreme_Lw() const { ASSERT(local_extreme_ != ExtremeType::unknown); return extreme_Lw_; }
  void set_opposite_direction_w(double opposite_direction_w) const { opposite_direction_w_ = opposite_direction_w; }
  double opposite_direction_w() const { ASSERT(local_extreme_ != ExtremeType::unknown); return opposite_direction_w_; }

  void explored(HorizontalDirection hdirection) const
  {
    ASSERT(hdirection != HorizontalDirection::undecided);
    int explore_flag = hdirection == HorizontalDirection::left ? explored_left : explored_right;
    explored_ |= explore_flag;
  }

  bool done() const
  {
    return explored_ == (explored_left|explored_right);
  }

  bool is_explored(HorizontalDirection hdirection) const
  {
    int explore_flag = hdirection == HorizontalDirection::left ? explored_left : explored_right;
    return (explored_ & explore_flag) != 0;
  }

  //---------------------------------------------------------------------------

#ifdef CWDEBUG
  double step() const
  {
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return next_sample_->w() - w();
  }

  void print_on(std::ostream& os) const
  {
    os << "{";
    Sample::print_on(os);
    if (explored_)
    {
      os << " E:";
      if (explored_ == explored_left)
        os << "←";
      else if (explored_ == explored_right)
        os << "→";
      else
        os << "↔";
    }
    os << ", cubic:";
    if (type_ == CubicToNextSampleType::unknown)
      os << "<none>";
    else
      os << cubic_ << ", type:" << type_;
    os << "}";
  }
#endif
};

} // namespace gradient_descent
