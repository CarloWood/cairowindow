#pragma once

#include "SampleNode.h"
#include <list>

namespace gradient_descent {

class ExtremeChain
{
 private:
  SampleNode::list_type sample_node_list_;      // Sorted list with all samples.
  SampleNode::iterator last_;                   // Points to the last sample that was added.

 public:
  void initialize(Sample&& first_sample);
  SampleNode::const_iterator insert(Sample&& new_sample);

  bool empty() const { return sample_node_list_.empty(); }
  SampleNode::const_iterator last() const { return last_; }
  SampleNode::const_iterator begin() const { return sample_node_list_.begin(); }
  SampleNode::const_iterator end() const { return sample_node_list_.end(); }
};

} // namespace gradient_descent
