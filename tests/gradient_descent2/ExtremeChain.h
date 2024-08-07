#pragma once

#include "SampleNode.h"
#include <list>

namespace gradient_descent {

class ExtremeChain
{
 public:
  static constexpr double negligible_scale_fraction = 1e-6;     // Smaller than this is ignored as non-existent.
  static constexpr double significant_scale_fraction = 1e-5;    // Assumed larger than float-point round-off errors.

 private:
  SampleNode::list_type sample_node_list_;      // Sorted list with all samples.
  SampleNode::const_iterator last_;             // Points to the last sample that was added.
  SampleNode::const_iterator larger_;           // Points to the first sample that is larger than the w value passed to find_larger (or end()).
#ifdef CWDEBUG
  double new_w_;                                // Copy of the value passed to the last call to find_larger.
#endif

 public:
  void initialize(Sample&& first_sample);
  bool find_larger(double new_w);
  SampleNode::const_iterator insert(Sample&& new_sample);
  void reuse(SampleNode::const_iterator const& sample_node) { last_ = sample_node; }

  bool empty() const { return sample_node_list_.empty(); }
  std::pair<SampleNode::const_iterator, bool> duplicate(double scale) const;
  SampleNode::const_iterator last() const { return last_; }
  SampleNode::const_iterator begin() const { return sample_node_list_.begin(); }
  SampleNode::const_iterator end() const { return sample_node_list_.end(); }
};

} // namespace gradient_descent
