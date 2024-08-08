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
  using list_type = std::list<SampleNode>;              // The actual list is member of ExtremeChain.
  using iterator = list_type::iterator;
  using const_iterator = list_type::const_iterator;

#if CW_DEBUG
 public:
  // Return value of find_extreme when no extreme is available.
  static constexpr double uninitialized_magic = 12345678.876543211;
#endif

 private:
  mutable math::CubicPolynomial cubic_;         // The cubic that fits this and the next Sample.
  mutable CubicToNextSampleType type_;          // The type of this cubic.
  mutable Scale scale_;                         // The scale that belongs to this cubic.
  mutable const_iterator next_sample_;          // The right-sample used for cubic_, copy of what was passed to initialize_cubic.

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
  bool has_extreme(ExtremeType extreme_type) const
  {
    ASSERT(extreme_type != ExtremeType::unknown);
    return extreme_type == ExtremeType::minimum ? has_minimum(type_) : has_maximum(type_);
  }

  // Accessors.
  math::CubicPolynomial const& cubic() const { return cubic_; }
  CubicToNextSampleType type() const { return type_; }
  Scale const& scale() const { return scale_; }
  bool is_rising() const
  {
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return (static_cast<int>(type_) & rising_bit) != 0 ||
      ((static_cast<int>(type_) & (minimum_bit|maximum_bit)) != 0 && next_sample_->Lw() > Lw());
  }
  bool is_falling() const
  {
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return (static_cast<int>(type_) & falling_bit) != 0 ||
      ((static_cast<int>(type_) & (minimum_bit|maximum_bit)) != 0 && next_sample_->Lw() < Lw());
  }
  double w_scale_estimate() const { ASSERT(next_sample_->w() > w()); return next_sample_->w() - w(); }
  double Lw_scale_estimate() const { return std::abs(next_sample_->Lw() - Lw()); }
  double dLdw_scale_estimate() const { return std::abs((next_sample_->Lw() - Lw()) / (next_sample_->w() - w())); }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{";
    Sample::print_on(os);
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
