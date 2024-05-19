#pragma once

#include "QuadraticPolynomial.h"
#include "Scale.h"
#include "Sample.h"
#include <array>

namespace gradient_descent {
using utils::has_print_on::operator<<;

enum HorizontalDirection
{
  left = -1,
  unknown_horizontal_direction = 0,
  right = 1
};

inline HorizontalDirection opposite(HorizontalDirection hdirection)
{
  return static_cast<HorizontalDirection>(-hdirection);
}

std::string to_string(HorizontalDirection hdirection);

inline std::ostream& operator<<(std::ostream& os, HorizontalDirection hdirection)
{
  return os << to_string(hdirection);
}

class Approximation
{
 private:
  int number_of_relevant_samples_{0};                   // The number of valid samples in relevant_samples_, used for this approximation.
  int current_index_{1};                                // The index of the last sample that was added (iff number_of_relevant_samples_ > 0).
  std::array<Sample const*, 2> relevant_samples_;       // Pointers to up to two samples that take part in this approximation.
  math::QuadraticPolynomial parabola_;                  // A linear or parabolic approximation.
  Scale parabola_scale_;                                // A measure of over what interval the parabolic approximation was tested to be correct.

#ifdef CWDEBUG
  bool is_part_of_extreme_{false};                      // Set when this is a LocalExtreme::approximation_.
#endif

 public:
  // Called with the latest samples that are expected to match this parabola (or that
  // should construct the parabola when there are not already two relevant samples stored).
  ScaleUpdate add(Sample const* current, bool update_scale_only, bool current_is_replacement);

#ifdef CWDEBUG
  void set_is_part_of_extreme()
  {
    is_part_of_extreme_ = true;
  }
#endif

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
    os << "{parabola:" << parabola_ << ", parabola_scale:" << parabola_scale_ << "}";
  }
#endif
};

} // namespace gradient_descent
