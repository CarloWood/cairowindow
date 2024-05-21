#pragma once

#include "Scale.h"
#include "Sample.h"
#include "utils/Array.h"
#include "debug.h"

namespace gradient_descent {

struct HistoryIndexCategory;
using HistoryIndex = utils::ArrayIndex<HistoryIndexCategory>;

class History
{
 public:
  static constexpr int size = 9;

 private:
  utils::Array<Sample, size, HistoryIndex> samples_;
  HistoryIndex current_;                // Index of the last Sample that was added.
  HistoryIndex prev_;                   // Index to the Sample that was added before current.
  int old_samples_ = 0;                 // Samples elsewhere that should not be used for fitting a polynomial.
  int total_number_of_samples_ = 0;     // The number of samples taken (not necessarily equal to the number that is (still) in the history).

 public:
  HistoryIndex add(double w, double Lw, double dLdw, Scale const& scale, bool& current_is_replacement)
  {
    DoutEntering(dc::notice, "History::add(" << w << ", " << Lw << ", " << dLdw << ", " << scale << ", current_is_replacement)");

    // If the new sample is very close to the current one, don't add it, just replace the current sample.
    bool almost_equal = !current_.undefined() && scale.negligible(samples_[current_].w() - w);
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
      current_ = HistoryIndex{(prev_ + 1).get_value() % size};
      ++total_number_of_samples_;
    }

    samples_[current_].set_values(w, Lw, dLdw);
    return current_;
  }

  void reset()
  {
    old_samples_ = total_number_of_samples_;
  }

  int relevant_samples() const
  {
    return total_number_of_samples_ - old_samples_;
  }

  static HistoryIndex before(HistoryIndex i)
  {
    ASSERT(!i.undefined() && !i.is_zero());
    return HistoryIndex{(i + size - 1).get_value() % size};
  }

  Sample const& operator[](HistoryIndex i) const
  {
    ASSERT(!i.undefined());
    ASSERT(i.get_value() < total_number_of_samples_);
    return samples_[i];
  }
  Sample const& current() const
  {
    ASSERT(!current_.undefined());
    ASSERT(current_.get_value() < total_number_of_samples_);
    return samples_[current_];
  }
  Sample const& prev(int count) const
  {
    ASSERT(total_number_of_samples_ > count);
    HistoryIndex index{(current_.get_value() + size - count) % size};
    return samples_[index];
  }

  int total_number_of_samples() const { return total_number_of_samples_; }
};

} // namespace gradient_descent
