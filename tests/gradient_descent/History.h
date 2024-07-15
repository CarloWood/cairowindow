#pragma once

#include "Scale.h"
#include "Sample.h"
#include "AlgorithmEventType.h"
#include "HistoryIndex.h"
#include "utils/Array.h"
#include <functional>
#include "debug.h"

namespace gradient_descent {

class History
{
 public:
  static constexpr int size = 9;

 private:
  utils::Array<Sample, size, HistoryIndex> samples_;
  HistoryIndex current_;                // Index of the last Sample that was added.
  HistoryIndex prev_;                   // Index to the Sample that was added before current.
  int old_samples_ = 0;                 // Samples elsewhere that should not be used for fitting a cubic.
  int total_number_of_samples_ = 0;     // The number of samples taken (not necessarily equal to the number that is (still) in the history).

#ifdef CWDEBUG
  events::Server<AlgorithmEventType>& event_server_;
#endif

 public:
#ifdef CWDEBUG
  History(events::Server<AlgorithmEventType>& event_server) : event_server_(event_server) { }
#endif

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

#ifdef CWDEBUG
    event_server_.trigger(AlgorithmEventType{history_add_event, current_, samples_[current_], std::to_string(total_number_of_samples_ - 1)});
#endif

    return current_;
  }

  // Return to, unless there are samples between from and to, in that case return the sample that is closest to from.
  // Samples for which `ignore` returns true must be ignored (we never clamp on those).
  double clamp(double from, double to, std::function<bool(Sample const&)> const& ignore, HistoryIndex& clamped_history_index)
  {
    clamped_history_index.set_to_undefined();
    HistoryIndex const iend = std::min(HistoryIndex{total_number_of_samples_}, samples_.iend());
    for (HistoryIndex index = samples_.ibegin(); index != iend; ++index)
    {
      if (ignore(samples_[index]))
        continue;
      double const sample_w = samples_[index].w();
      if ((from < sample_w && sample_w < to) || (to < sample_w && sample_w < from))
      {
        to = sample_w;
        clamped_history_index = index;
      }
    }
    return to;
  }

  void reset()
  {
    DoutEntering(dc::notice, "History::reset()");
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
