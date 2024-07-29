#pragma once

#include "Sample.h"
#include "CubicToNextSampleType.h"
#include "HorizontalDirection.h"
#include "ExtremeType.h"
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
#if CW_DEBUG
 public:
  // Return value of find_extreme when no extreme is available.
  static constexpr double uninitialized_magic = 12345678.876543211;
#endif

 private:
  mutable math::CubicPolynomial cubic_;         // The cubic that fits this and the next Sample.
  mutable CubicToNextSampleType type_;          // The type of this cubic.

 public:
  SampleNode(Sample&& sample) : Sample(std::move(sample)), type_(CubicToNextSampleType::unknown) { }

  void initialize_cubic(SampleNode const& next
      COMMA_CWDEBUG_ONLY(events::Server<AlgorithmEventType>& event_server, bool this_is_last)) const;

  double find_extreme(Sample const& next, Region& region_out, ExtremeType& extreme_type) const;

  // Accessors.
  math::CubicPolynomial const& cubic() const { return cubic_; }
  CubicToNextSampleType type() const { return type_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{";
    Sample::print_on(os);
    os << ", cubic:" << cubic_ << ", type:" << type_ << "}";
  }
#endif
};

} // namespace gradient_descent
