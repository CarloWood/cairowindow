#pragma once

#include "Scale.h"
#include <array>
#include "debug.h"

namespace gradient_descent {

template<ConceptSample T>
class History
{
 public:
  static constexpr int size = 9;

 private:
  std::array<T, size> samples_;
  int current_ = -1;                    // Index of the last Sample that was added.
  int prev_ = -1;                       // Index to the Sample that was added before current.
  int old_samples_ = 0;                 // Samples elsewhere that should not be used for fitting a polynomial.
  int total_number_of_samples_ = 0;     // The number of samples taken (not necessarily equal to the number that is (still) in the history).

 public:
  T const& add(double w, double Lw, double dLdw, Scale const& scale, bool& current_is_replacement)
  {
    // If the new sample is very close to the current one, don't add it, just replace the current sample.
    bool almost_equal = (current_ != -1 && std::abs(samples_[current_].w() - w) < 0.001 * scale.or_zero());
    if (almost_equal)
    {
      Dout(dc::notice, "Replacing history point " << (total_number_of_samples_ - 1) << " at " << samples_[current_].w() << " --> " << w);
      current_is_replacement = true;
    }
    else
    {
      Dout(dc::notice, "Appending to history: point " << total_number_of_samples_ << " at w = " << w);
      current_is_replacement = false;

      prev_ = current_;
      current_ = (prev_ + 1) % size;
      ++total_number_of_samples_;
    }

    samples_[current_].set_values(w, Lw, dLdw);
    return samples_[current_];
  }

  void reset()
  {
    old_samples_ = total_number_of_samples_;
  }

  int relevant_samples() const
  {
    return total_number_of_samples_ - old_samples_;
  }

  static int before(int i) { ASSERT(0 <= i); return (i + size - 1) % size; }

  T const& current() const { ASSERT(current_ != -1); ASSERT(current_ < total_number_of_samples_); return samples_[current_]; }
  T const& prev() const { ASSERT(prev_ != -1); ASSERT(prev_ < total_number_of_samples_); return samples_[prev_]; }
  T const& prev_prev() const { ASSERT(prev_ != -1); ASSERT(prev_ < total_number_of_samples_); return samples_[before(prev_)]; }

  int total_number_of_samples() const { return total_number_of_samples_; }
};

} // namespace gradient_descent
