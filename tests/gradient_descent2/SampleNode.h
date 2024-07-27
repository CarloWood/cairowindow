#pragma once

#include "Sample.h"
#include "CubicToNextSampleType.h"
#include <memory>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include <iostream>
#endif

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

// SampleNode's form a linked list, stored in Algorithm, of all Sample's done
// so far. They are strictly ordered by Sample::w_ value: the w value of prev_
// is less than, and the w value of next_ is larger than that of this sample.
//
// Any additional data stored in SampleNode refers to the cubic that fits
// through this Sample and the next_ Sample.
class SampleNode : public Sample
{
 private:
  mutable CubicToNextSampleType type_;          // The type of cubic that fits this and the next Sample.

 public:
  SampleNode(Sample&& sample) : Sample(std::move(sample)), type_(CubicToNextSampleType::unknown) { }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{";
    Sample::print_on(os);
    os << ", type:" << type_ << "}";
  }
#endif
};

} // namespace gradient_descent
