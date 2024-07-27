#pragma once

#include "SampleNode.h"
#include <list>

namespace gradient_descent {

class ExtremeChain
{
 public:
  using sample_node_list_type = std::list<SampleNode>;
  using iterator = sample_node_list_type::iterator;
  using const_iterator = sample_node_list_type::const_iterator;

 private:
  sample_node_list_type sample_node_list_;      // Sorted list with all samples.
  iterator last_;                               // Points to the last sample that was added.

 public:
  void initialize(Sample&& first_sample);
  ExtremeChain::const_iterator insert(Sample&& new_sample);

  bool empty() const { return sample_node_list_.empty(); }
  SampleNode const& last() const { return *last_; }
  ExtremeChain::const_iterator begin() const { return sample_node_list_.begin(); }
  ExtremeChain::const_iterator end() const { return sample_node_list_.end(); }
};

} // namespace gradient_descent
