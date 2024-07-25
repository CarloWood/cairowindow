#pragma once

#include "Scale.h"
#include "Sample.h"
#include "HorizontalDirection.h"
#include "ExtremeType.h"
#include "../QuadraticPolynomial.h"
#include "../CubicPolynomial.h"
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
  std::array<Sample const*, 2> relevant_samples_;       // Pointers to up to two samples that take part in this approximation. If
                                                        // number_of_relevant_samples_ is two then one with the smallest w is stored at index 0.
  Sample back_tracking_pivot_;                          // Copy of a sample that is used while back tracking and might be needed longer than
                                                        // the lifetime of the original in the history.
  math::CubicPolynomial cubic_;                         // A cubic approximation, only valid if number_of_relevant_samples_ == 2.
  Scale scale_;                                         // A measure of over what interval the approximation was tested to be correct.
  bool is_extreme_{false};                              // Set when this is a LocalExtreme::approximation_.
  bool already_had_two_relevant_samples_{false};        // Used as carry between add and update_scale.

 private:
  // Only called by the move constructor.
  void reset()
  {
    // Forget any previous samples.
    number_of_relevant_samples_ = 0;
    current_index_ = 1;
    scale_.reset();
  }

 public:
  Approximation() = default;
  Approximation(Approximation const&) = default;

  Approximation(Approximation&& orig) : number_of_relevant_samples_{orig.number_of_relevant_samples_}, current_index_{orig.current_index_},
    relevant_samples_{orig.relevant_samples_}, scale_{std::move(orig.scale_)},
    is_extreme_{orig.is_extreme_}, cubic_{std::move(orig.cubic_)}, already_had_two_relevant_samples_{orig.already_had_two_relevant_samples_}
  {
    orig.reset();
  }

  Approximation& operator=(Approximation const&) = default;
  Approximation& operator=(Approximation&& orig)
  {
    number_of_relevant_samples_ = orig.number_of_relevant_samples_;
    current_index_ = orig.current_index_;
    relevant_samples_ = orig.relevant_samples_;
    scale_ = std::move(orig.scale_);
    is_extreme_ = orig.is_extreme_;
    cubic_ = std::move(orig.cubic_);
    already_had_two_relevant_samples_ = orig.already_had_two_relevant_samples_;
    orig.reset();
    return *this;
  }

  // Called with the latest samples that are expected to match this cubic approximation (or that
  // should construct the cubic when there are not already two relevant samples stored).
  void add(Sample const* current, bool current_is_replacement, ExtremeType next_extreme_type, bool back_tracking);
  ScaleUpdate update_scale(bool current_is_replacement, ExtremeType next_extreme_type);
  ScaleUpdate update_local_extreme_scale(Sample const& current);

  Weight find_extreme(Region& region, ExtremeType& extreme_type) const;
  void set_current_index(Region region);

  Sample const& current() const { ASSERT(number_of_relevant_samples_ > 0); return *relevant_samples_[current_index_]; }
  Sample const& prev() const { ASSERT(number_of_relevant_samples_ > 1); return *relevant_samples_[1 - current_index_]; }

  void set_is_extreme()
  {
    is_extreme_ = true;
  }

  void unset_is_extreme()
  {
    is_extreme_ = false;
  }

  int number_of_relevant_samples() const { return number_of_relevant_samples_; }
  Scale const& scale() const { return scale_; }
  bool is_extreme() const { return is_extreme_; }
  double at(double w) const { return cubic_(w); }
  math::CubicPolynomial const& cubic() const { return cubic_; }
  HorizontalDirection prev_to_current() const { return prev().w() > current().w() ? HorizontalDirection::left : HorizontalDirection::right; }

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
    os << "cubic:" << cubic_ <<
        ", scale:" << scale_;
    if (number_of_relevant_samples_ > 1)
    {
      std::array<double, 2> extremes;
      Debug(libcw_do.off());
      int n = cubic_.get_extremes(extremes);
      Debug(libcw_do.on());
      if (n == 0)
        os << " [no extremes]";
      else if (n == 1)
        os << " [derivative root = " << extremes[0] << "]";
      else
      {
        int index_minimum = (cubic_[3] > 0.0) ? 1 : 0;
        os << " [x_minimum = " << extremes[index_minimum] << ", x_maximum = " << extremes[1 - index_minimum] << "]";
      }
    }
    os << "}";
  }
#endif
};

} // namespace gradient_descent
