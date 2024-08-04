#include "sys.h"
#include "ExtremeChain.h"

namespace gradient_descent {

void ExtremeChain::initialize(Sample&& first_sample)
{
  last_ = sample_node_list_.emplace(sample_node_list_.begin(), std::move(first_sample));
}

SampleNode::const_iterator ExtremeChain::insert(Sample&& new_sample)
{
  // Use initialize() for the first sample.
  ASSERT(!empty());

  double const new_w = new_sample.w();

  // Find the two adjacent samples that lay left and right of new_sample.

  // Start at last_, because it is very likely that new_w needs to be inserted
  // right next to it, or at least in the neighborhood.

  // This is going to point to the first element that is larger than new_w (or end() if none such exists).
  SampleNode::const_iterator larger = last_;

  if (larger->w() > new_w)
  {
    // If the initial value is larger then step to the left until that is no longer larger.
    while (larger != sample_node_list_.begin())
    {
      if ((--larger)->w() <= new_w)     // Not larger anymore?
      {
        ++larger;                       // Undo the decrement and exit.
        break;
      }
    }
  }
  else
  {
    // If the initial value is smaller then step to the right until we find a sample that is larger.
    while (++larger != sample_node_list_.end() && larger->w() <= new_w)
      ;
  }

  // Sanity check for the above code.
  ASSERT(larger == sample_node_list_.end() || larger->w() > new_w);
  ASSERT(larger == sample_node_list_.begin() || std::prev(larger)->w() <= new_w);

  // Insert the new sample before `larger`.
  last_ = sample_node_list_.insert(larger, std::move(new_sample));
  return last_;
}

} // namespace gradient_descent
