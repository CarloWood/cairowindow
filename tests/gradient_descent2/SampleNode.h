#pragma once

#include "LocalExtreme.h"
#include "Sample.h"
#include "Scale.h"
#include "CubicToNextSampleType.h"
#include "HorizontalDirection.h"
#include "ExtremeType.h"
#include "CriticalPointType.h"
#include "math/CubicPolynomial.h"
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
  const_iterator next_node_;                            // The right-sample, if any, used for cubic_,
                                                        // copy of what was passed to initialize_cubic.
  math::CubicPolynomial cubic_;                         // The cubic that fits this and next_node_, if the latter isn't chain_.end().
#if CW_DEBUG
  bool cubic_initialized_{false};                       // Set to true after cubic_ was initialized.
#endif
  CubicToNextSampleType type_;                          // The type of cubic_.
  mutable Scale scale_;                                 // The scale that belongs to this cubic.

  mutable std::unique_ptr<LocalExtreme> local_extreme_; // Points to LocalExtreme instance iff this cubic was used for a local extreme.

 public:
  SampleNode(Sample&& sample) : Sample(std::move(sample)), type_(CubicToNextSampleType::unknown) { }

  void initialize_cubic(const_iterator next
      COMMA_CWDEBUG_ONLY(ExtremeType next_extreme_type, events::Server<AlgorithmEventType>& event_server, bool this_is_last));

  void change_type_to_left_extreme(ExtremeType extreme_type);
  void change_type_to_right_extreme(ExtremeType extreme_type);

 public:
  double find_extreme(Sample const& next, ExtremeType& extreme_type) const;

  const_iterator next_node() const
  {
    // Call initialize_cubic before using this member function.
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return next_node_;
  }

  void set_scale(CriticalPointType type, double critical_point_w, double left_edge_w, double right_edge_w) const
  {
    scale_.set(type, critical_point_w, left_edge_w, right_edge_w, cubic_.inflection_point());
  }

  //---------------------------------------------------------------------------
  // Accessors for the cubic.

  // Returns the cubic that was fitted between this and the next sample.
  math::CubicPolynomial const& cubic() const { ASSERT(cubic_initialized_); return cubic_; }

  // Returns the type of the cubic that was fitted between this and the next sample.
  CubicToNextSampleType type() const { return type_; }

  // Returns true if there is a node on the right and the cubic in between is initialized.
  bool is_cubic() const { return type_ != CubicToNextSampleType::unknown; }

  // The scale of this cubic.
  Scale const& scale() const { return scale_; }

  // Returns true if this type falls in the 'rising' class.
  bool is_rising() const
  {
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return (static_cast<int>(type_) & rising_bit) != 0 ||
      ((static_cast<int>(type_) & (minimum_bit|maximum_bit)) != 0 && next_node_->Lw() > Lw());
  }

  // Returns true if this type falls in the 'falling' class.
  bool is_falling() const
  {
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return (static_cast<int>(type_) & falling_bit) != 0 ||
      ((static_cast<int>(type_) & (minimum_bit|maximum_bit)) != 0 && next_node_->Lw() < Lw());
  }

  // Returns true if the cubic fitted between this and the next sample has the given extreme in between the two samples.
  bool has_unfound_extreme(ExtremeType extreme_type) const
  {
    ASSERT(extreme_type != ExtremeType::unknown);
    return extreme_type == ExtremeType::minimum ? has_unfound_minimum(type_) : has_unfound_maximum(type_);
  }

  bool has_extreme(ExtremeType extreme_type) const
  {
    ASSERT(extreme_type != ExtremeType::unknown);
    return extreme_type == ExtremeType::minimum ? has_minimum(type_) : has_maximum(type_);
  }

  bool has_extreme() const
  {
    return (static_cast<int>(type_) & (minimum_bit|left_min_bit|right_min_bit|maximum_bit|left_max_bit|right_max_bit)) != 0;
  }

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

  bool has_flat_extreme(ExtremeType extreme_type, SideOfCubic on_side) const
  {
    int bits_to_test = left_min_bit << (static_cast<int>(extreme_type) + 1 + on_side);
    return (static_cast<int>(type_) & bits_to_test) != 0;
  }

  //---------------------------------------------------------------------------
  // Local extreme data.

  // Returns the w value of the underlying extreme.
  double extreme_w() const
  {
    // Call initialize_cubic before using this member function.
    ASSERT(type_ != CubicToNextSampleType::unknown);
    // Call set_local_extreme before calling this function.
    ASSERT(local_extreme_);
    return scale_.critical_point_w();
  }

  //---------------------------------------------------------------------------
  // Scale estimates that can be used when scale is not available yet.
  double w_scale_estimate() const
  {
    ASSERT(next_node_->w() > w());
    return next_node_->w() - w();
  }
  double Lw_scale_estimate() const
  {
    ASSERT(cubic_initialized_);
    // The derivaties of the cubic in the current point (w()).
    double const dQ = dLdw();
    double const ddQ = cubic_.second_derivative(w());
    double const dddQ = 6.0 * cubic_[3];
    double const u = next_node_->w() - w();

    // The RMS of Q'(x) over the interval between this and the next sample. See rms.m.
    double const avg_absolute_derivative =
      std::sqrt(utils::square(dQ) + (dQ * ddQ +
            ((dQ * dddQ + utils::square(ddQ)) / 3.0 + (0.25 * ddQ * dddQ + 0.05 * utils::square(dddQ) * u) * u) * u) * u);

    return avg_absolute_derivative;
  }
  double dLdw_scale_estimate() const { return Lw_scale_estimate() / w_scale_estimate(); }

  //---------------------------------------------------------------------------
  // Local extreme functions.

  // Mark this node as the left node of a cubic whose local_extreme type is a local extreme of the underlying function L(w).
  void set_local_extreme(ExtremeType local_extreme, const_iterator chain_end) const
  {
    // Pass minimum or maximum.
    ASSERT(local_extreme != ExtremeType::unknown);
    // Call set_scale before calling this function.
    ASSERT(scale_.type() != CriticalPointType::none);
    // Only call this function once.
    ASSERT(!local_extreme_);
    local_extreme_ = std::make_unique<LocalExtreme>(local_extreme, cubic_(scale_.critical_point_w()), chain_end
        COMMA_CWDEBUG_ONLY('[' + std::to_string(label()) + ']'));
  }

  // Return true iff set_local_extreme was called.
  bool is_local_extreme() const { return static_cast<bool>(local_extreme_); }

  LocalExtreme& local_extreme() const
  {
    // Only call this function if is_local_extreme returns true.
    ASSERT(local_extreme_);
    return *local_extreme_;
  }

  //---------------------------------------------------------------------------

#ifdef CWDEBUG
  double step() const
  {
    ASSERT(type_ != CubicToNextSampleType::unknown);
    return next_node_->w() - w();
  }

  void print_on(std::ostream& os) const
  {
    os << "{";
    Sample::print_on(os);
    if (local_extreme_)
      local_extreme_->print_on(os);
    os << ", cubic:";
    if (!cubic_initialized_)
      os << "<unknown>";
    else
      os << cubic_ << ", type:" << type_;
    os << "}";
  }
#endif
};

} // namespace gradient_descent
