#include "sys.h"
#include "ExtremeChain.h"
#ifdef CWDEBUG
#include "Algorithm.h"
#endif

namespace gradient_descent {

void ExtremeChain::initialize(Sample&& first_sample)
{
  last_ = sample_node_list_.emplace(sample_node_list_.begin(), std::move(first_sample));
}

void ExtremeChain::find_larger(double const new_w)
{
  DoutEntering(dc::notice, "ExtremeChain::find_larger(" << new_w << ")");

  // Use initialize() for the first sample.
  ASSERT(!empty());

  // Find the two adjacent samples that lay left and right of new_sample.

  // Start at last_, because it is very likely that new_w needs to be inserted
  // right next to it, or at least in the neighborhood.

  // This is going to point to the first element that is larger than new_w (or end() if none such exists).
  larger_ = last_;

  if (larger_->w() > new_w)
  {
    // If the initial value is larger then step to the left until that is no longer larger.
    while (larger_ != sample_node_list_.begin())
    {
      if ((--larger_)->w() <= new_w)    // Not larger anymore?
      {
        ++larger_;                      // Undo the decrement and exit.
        break;
      }
    }
  }
  else
  {
    // If the initial value is smaller (or equal) then step to the right until we find a sample that is larger.
    while (++larger_ != sample_node_list_.end() && larger_->w() <= new_w)
      ;
  }

  // Sanity check for the above code.
  ASSERT(larger_ == sample_node_list_.end() || larger_->w() > new_w);
  ASSERT(larger_ == sample_node_list_.begin() || std::prev(larger_)->w() <= new_w);
  Debug(new_w_ = new_w);

#ifdef CWDEBUG
  if (larger_ == sample_node_list_.end())
    Dout(dc::notice, "None of the samples in the list are larger than " << new_w);
  else
    Dout(dc::notice, "Found larger = " << *larger_);
#endif
}

std::pair<SampleNode::const_iterator, bool> ExtremeChain::duplicate(double scale) const
{
  std::pair<SampleNode::const_iterator, bool> ibp{{}, false};
  double const significant_difference = significant_scale_fraction * scale;
  if (AI_UNLIKELY(larger_ != sample_node_list_.begin() && std::abs(std::prev(larger_)->w() - new_w_) < significant_difference))
    ibp = {std::prev(larger_), true};
  else if (AI_UNLIKELY(larger_ != sample_node_list_.end() && std::abs(larger_->w() - new_w_) < significant_difference))
    ibp = {larger_, true};
  return ibp;
}

SampleNode::iterator ExtremeChain::insert(Sample&& new_sample)
{
  // Call find_larger with the w value of the next sample.
  ASSERT(new_w_ == new_sample.w());

  // Insert the new sample before `larger_`.
  auto new_node = sample_node_list_.insert(larger_, std::move(new_sample));
  last_ = new_node;
  return new_node;
}

#ifdef CWDEBUG
void ExtremeChain::dump(Algorithm const* algorithm) const
{
  Dout(dc::notice, "Current chain:");
  NAMESPACE_DEBUG::Indent indent(2);

  SampleNode::const_iterator node = sample_node_list_.begin();
  SampleNode::const_iterator end = sample_node_list_.end();

  // Print at most 7 nodes around the last node that was added.
  if (sample_node_list_.size() > 7)
  {
    end = node = algorithm->debug_chain().last();
    int left_count = 0;
    while (left_count < 3 && node != sample_node_list_.begin())
    {
      --node;
      ++left_count;
    }
    int right_count = 0;
    while (right_count < 7 - left_count && end != sample_node_list_.end())
    {
      ++end;
      ++right_count;
    }
    while (left_count + right_count < 7 && node != sample_node_list_.begin())
    {
      --node;
      ++left_count;
    }
  }

  while (node != end)
  {
    Dout(dc::notice|continued_cf, '[' << node->label() << "] " << std::setprecision(12) << node->w());
    if (node->is_local_extreme())
      Dout(dc::continued, " [" << node->get_extreme_type() << ']');
    Dout(dc::finish, " " << node->scale());
    if (node->type() != CubicToNextSampleType::unknown)
    {
      NAMESPACE_DEBUG::Indent indent(2);
      Dout(dc::notice|continued_cf, std::left << std::setw(5) << to_utf8_art(node->type()) << '+' << std::setw(10) << node->step() << node->cubic());
      if (node == algorithm->debug_cubic_used())
        Dout(dc::continued, " [cubic_used_]");
      if (node == algorithm->debug_left_of())
        Dout(dc::continued, " [left_of]");
      if (node == algorithm->debug_right_of())
        Dout(dc::continued, " [right_of]");
      Dout(dc::finish, "");
    }
    ++node;
  }
}

namespace {

enum CubicEndShape
{
  flat_high,
  flat_low,
  downhill,
  uphill,
  plus_inf
};

CubicEndShape get_end(CubicToNextSampleType type, bool left)
{
  using enum CubicToNextSampleType;
  switch (type)
  {
    case flat:                  // __
      return flat_low;
    case up:                    // /
      return uphill;
    case down:                  // \.
      return downhill;
    case right_stop:            // /^
      return left ? uphill : plus_inf;
    case left_stop:             // ^\.
      return left ? plus_inf : downhill;
    case right_min:             // _/
      return left ? flat_low : uphill;
    case left_min:              // \_
      return left ? downhill : flat_low;
    case right_max:             // ‾\.
      return left ? flat_high : downhill;
    case left_max:              // /‾
      return left ? uphill : flat_high;
    case right_max_left_min:    // ‾\_
      return left ? flat_high : flat_low;
    case right_min_left_max:    // _/‾
      return left ? flat_low : flat_high;
    case min:                   // \/
      return left ? downhill : uphill;
    case right_max_min:         // ‾\/
      return left ? flat_high : uphill;
    case min_left_max:          // \/‾
      return left ? downhill : flat_high;
    case min_max:               // \/\.
      return downhill;
    case max_min:               // /\/
      return uphill;
    case max:                   // /\.
      return left ? uphill : downhill;
    case max_left_min:          // /\_
      return left ? uphill : flat_low;
    case right_min_max:         // _/\.
      return left ? flat_low : downhill;
    default:
      ASSERT(false);
  }
  AI_NEVER_REACHED
}

} // namespace

void ExtremeChain::sanity_check(Algorithm const* algorithm) const
{
  if (sample_node_list_.empty())
    return;
  CubicToNextSampleType prev_type = CubicToNextSampleType::unknown;
  SampleNode::const_iterator last_extreme;
  ExtremeType last_extreme_type = ExtremeType::unknown;     // Extreme type of the last extreme;
  // Run over all node pairs {node, next} that make up cubics.
  SampleNode::const_iterator node = sample_node_list_.begin();
  for (SampleNode::const_iterator next = std::next(node); next != sample_node_list_.end(); prev_type = node->type(), node = next++)
  {
    // The type of the cubic between every two samples in the chain must be known.
    ASSERT(node->type() != CubicToNextSampleType::unknown);

    bool is_first = prev_type == CubicToNextSampleType::unknown;
    bool is_last = std::next(next) == sample_node_list_.end();

    // Only the very first one may be a left_stop.
    ASSERT(!is_first || node->type() != CubicToNextSampleType::left_stop);
    // Only the last one may be a right_stop.
    ASSERT(!is_last || node->type() != CubicToNextSampleType::right_stop);

    // Does this cubic contain an extreme?
    if (node->is_local_extreme())
    {
      // Can't have two extrema of the same type on a row.
      if (last_extreme_type == node->get_extreme_type())
      {
        // This should never happen and probably means that we're trying to detect a local extreme that was already found before.
        DoutFatal(dc::core, "The node at " << static_cast<Sample const&>(*last_extreme) << " has the same type as the node at " <<
            static_cast<Sample const&>(*node) << " (" << last_extreme_type << ")!");
      }
      last_extreme_type = node->get_extreme_type();
      last_extreme = node;
    }

    // If this is the first cubic then we're done.
    if (is_first)
      continue;

    // The types must match.
    ASSERT(get_end(prev_type, false) == get_end(node->type(), true));
  }
}
#endif

} // namespace gradient_descent
