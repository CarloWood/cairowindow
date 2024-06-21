#pragma once

#include "Scale.h"
#include "Sample.h"
#include "../QuadraticPolynomial.h"
#include <array>

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Approximation
{
 private:
  int number_of_relevant_samples_{0};                   // The number of valid samples in relevant_samples_, used for this approximation.
  int current_index_{1};                                // The index of the last sample that was added (iff number_of_relevant_samples_ > 0).
  std::array<Sample const*, 2> relevant_samples_;       // Pointers to up to two samples that take part in this approximation.
  math::QuadraticPolynomial parabola_;                  // A linear or parabolic approximation.
  Scale parabola_scale_;                                // A measure of over what interval the parabolic approximation was tested to be correct.
  bool is_extreme_{false};                              // Set when this is a LocalExtreme::approximation_.

 public:
  // Called with the latest samples that are expected to match this parabola (or that
  // should construct the parabola when there are not already two relevant samples stored).
  ScaleUpdate add(Sample const* current, bool update_scale_only, bool current_is_replacement);

  Sample const& current() const { ASSERT(number_of_relevant_samples_ > 0); return *relevant_samples_[current_index_]; }
  Sample const& prev() const { ASSERT(number_of_relevant_samples_ > 1); return *relevant_samples_[1 - current_index_]; }

  void set_is_extreme()
  {
    is_extreme_ = true;
  }

  bool has_maximum() const
  {
    return parabola_[2] < 0.0;
  }

  void reset()
  {
    // Forget any previous samples.
    number_of_relevant_samples_ = 0;
    current_index_ = 1;
    parabola_scale_.reset();
  }

  int number_of_relevant_samples() const { return number_of_relevant_samples_; }
  math::QuadraticPolynomial const& parabola() const { return parabola_; }
  Scale const& parabola_scale() const { return parabola_scale_; }
  bool is_extreme() const { return is_extreme_; }

  int current_index() const
  {
    ASSERT(number_of_relevant_samples_ > 0);
    return current_index_;
  }

#ifdef CWDEBUG
  void print_based_on(std::ostream& os) const
  {
    os << "based on " << *relevant_samples_[current_index_];
    if (number_of_relevant_samples_ > 1)
      os << " and " << *relevant_samples_[1 - current_index_];
  }

  void print_on(std::ostream& os) const
  {
    os << "{parabola:" << parabola_;
    if (number_of_relevant_samples_ > 1)
      os << " [v_x = " << parabola_.vertex_x() << "]";
    os << ", parabola_scale:" << parabola_scale_ << "}";
  }
#endif
};

} // namespace gradient_descent
