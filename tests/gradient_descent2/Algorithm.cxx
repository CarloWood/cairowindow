#include "sys.h"
#include "SampleNode.h"
#include "Algorithm.h"
#include "IterationState.h"

namespace gradient_descent {

bool Algorithm::operator()(double& w, double Lw, double dLdw)
{
  DoutEntering(dc::notice, "Algorithm::operator()(" << w << ", " << Lw << ", " << dLdw << ")");

#ifdef CWDEBUG
  // Erase all previous curves (if they exist).
  event_server_.trigger(AlgorithmEventType{reset_event});
  if (have_expected_Lw_)
  {
    // Draw an arrow from where expected to end up, based on the approximation, and the actual value of L(w).
    event_server_.trigger(AlgorithmEventType{difference_event, w, expected_Lw_, Lw});
    have_expected_Lw_ = false;
  }
#endif

  // Create a new Sample for this sample.
  Sample current{w, Lw, dLdw COMMA_CWDEBUG_ONLY(label_context_)};
  Dout(dc::notice, "current = " << current);

#ifdef CWDEBUG
  // Draw the new sample.
  event_server_.trigger(AlgorithmEventType{new_sample_event, current});
#endif

  // Update kinetic energy. Returns false if too much energy was used.
  if (!update_energy(Lw))
    return handle_abort_hdirection(w);

  Dout(dc::notice, "Running " << state_);
  switch (state_)
  {
    case IterationState::initialization:
    {
      // We should only get here only the first time.
      ASSERT(chain_.empty());

      // Add the first sample.
      chain_.initialize(std::move(current));

#ifdef CWDEBUG
      {
        // Draw a line through this sample.
        math::QuadraticPolynomial line{Lw - dLdw * w, dLdw, 0.0};
        event_server_.trigger(AlgorithmEventType{quadratic_polynomial_event, line});
      }
#endif

      // Set w to the next value to probe. This uses chain_.last().
      handle_single_sample(w);

      state_ = IterationState::first_cubic;
      break;
    }

    case IterationState::first_cubic:
    {
      // Find the point where to insert the new sample. We do this by looking for the
      // first SampleNode (target) that has a w value that is greater than that of the
      // new sample and then inserting current in front of that.
      auto left_node = chain_.insert(std::move(current));
      auto right_node = std::next(left_node);
      // Fix the order.
      bool current_is_left = right_node != chain_.end();
      if (!current_is_left)
      {
        // In this case the iterator returned by insert is actually the right_node of the two.
        right_node = left_node;
        left_node = std::prev(right_node);
      }

      // The cubic is part of the left_node.
      left_node->initialize_cubic(*right_node COMMA_CWDEBUG_ONLY(event_server_, current_is_left));

      Region region;
      w = left_node->find_extreme(*right_node, region, next_extreme_type_);
      // w should be set if next_extreme_type_ was set.
      ASSERT(next_extreme_type_ == ExtremeType::unknown || w != SampleNode::uninitialized_magic);

      if (next_extreme_type_ == ExtremeType::unknown)
      {
        // We're going to look for a minimum.
        if (left_node->dLdw() > 0.0)
          left_of_ = &*left_node;
        else
          right_of_ = &*right_node;
      }
      else
      {
        if (w < left_node->w())
          left_of_ = &*left_node;
        else if (w <= right_node->w())
        {
          right_of_ = &*left_node;
          left_of_ = &*right_node;
        }
        else
          right_of_ = &*right_node;
      }
#ifdef CWDEBUG
      event_server_.trigger(AlgorithmEventType{left_of_right_of_event, left_of_, right_of_,
          next_extreme_type_, w, left_node->cubic()(w)});
#endif

      // Initialize the scale.
      {
        CriticalPointType scale_cp_type =
          next_extreme_type_ == ExtremeType::unknown ? CriticalPointType::inflection_point
                                                     : next_extreme_type_ == ExtremeType::minimum ? CriticalPointType::minimum
                                                                                                  : CriticalPointType::maximum;
        double critical_point_w =
          next_extreme_type_ == ExtremeType::unknown ? left_node->cubic().inflection_point()
                                                     : w;

        left_node->set_scale(scale_cp_type, critical_point_w, left_node, right_node);
#ifdef CWDEBUG
        event_server_.trigger(AlgorithmEventType{scale_draw_event, ScaleUpdate::initialized, *left_node, {}});
#endif
      }

      // If next_extreme_type_ was not set then this cubic has no extremes. In that case go for a minimum.
      if (next_extreme_type_ == ExtremeType::unknown)
        next_extreme_type_ = ExtremeType::minimum;

      state_ = IterationState::find_extreme;
      break;
    }

    case IterationState::find_extreme:
    {
      auto new_node = chain_.insert(std::move(current));
      auto right_node = std::next(new_node);

      if (right_node != chain_.end())
      {
        new_node->initialize_cubic(*right_node COMMA_CWDEBUG_ONLY(event_server_, true));        // true: called on new_node.
      }
      if (new_node != chain_.begin())
      {
        auto left_node = std::prev(new_node);
        left_node->initialize_cubic(*new_node COMMA_CWDEBUG_ONLY(event_server_, false));        // false: new_node is passed.
      }

      state_ = IterationState::next_sample;
      break;
    }

    case IterationState::next_sample:
    {
      w += bogus_;
      bogus_ = std::copysign(std::abs(bogus_) - 10.0, -bogus_);
      expected_Lw_ = 1800.0;
      have_expected_Lw_ = true;
      Debug(set_algorithm_str(w, "bogus"));
      break;
    }

    default:
      // Implement this state_.
      ASSERT(false);
  }

  // Ask for a new sample at `w` if not already finished.
  Dout(dc::notice, "Next probe: " << w);
  return state_ != IterationState::success;
}

void Algorithm::handle_single_sample(double& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_single_sample(" << w << ")");

  Sample const& sample = chain_.last();
  double step;
#ifdef CWDEBUG
  char const* algorithm_str;
#endif

  if (small_step_ == 0.0)       // Not defined yet?
  {
    step = -learning_rate_ * sample.dLdw();

    // Did we drop on a (local) extreme as a starting point?!
    if (std::abs(step) < epsilon || std::abs(step) < 1e-6 * std::abs(w))
    {
      if (hdirection_ == HorizontalDirection::undecided)
      {
        // In this case we can't do anything else than just make a step in some random direction.
        step = learning_rate_;
#ifdef CWDEBUG
        algorithm_str = "one sample, derivative is zero, hdirection is unknown";
#endif
      }
      else
      {
        step = static_cast<int>(hdirection_) * learning_rate_;
#ifdef CWDEBUG
        algorithm_str = "one sample, derivative is zero";
#endif
      }
    }
    else if (hdirection_ == HorizontalDirection::undecided)
    {
      // Just gradient descent: move downhill.
#ifdef CWDEBUG
      algorithm_str = "one sample, gradient descent";
#endif
    }
    else
    {
      // Make a step in the same horizontal direction.
      step = static_cast<int>(hdirection_) * std::abs(step);
#ifdef CWDEBUG
      algorithm_str = "one sample, same direction";
#endif
    }
  }
  else
  {
    // small_step_ should only be set once next_extreme_type_ has been set.
    ASSERT(next_extreme_type_ != ExtremeType::unknown);
    step = ((next_extreme_type_ == ExtremeType::minimum) == (sample.dLdw() > 0.0) ? -small_step_ : small_step_);
#ifdef CWDEBUG
    algorithm_str = "small step";
#endif
  }

#if 0 // FIXME: add scale back
  // This step could still be too small.
  if (...->scale().negligible(step))
  {
    // This wouldn't be good because then the new sample will replace
    // the current one and we'd still have just one sample.
    step = ...->scale().make_significant(step);
#ifdef CWDEBUG
    algorithm_str = "avoiding replacement";
#endif
  }
#endif
  w += step;
  have_expected_Lw_ = false;
  Debug(set_algorithm_str(w, algorithm_str));
}

#ifdef CWDEBUG
void Algorithm::set_algorithm_str(double new_w, char const* algorithm_str)
{
  algorithm_str_ = algorithm_str;
  Dout(dc::notice|continued_cf, std::setprecision(12) << chain_.last().w() << " --> " << chain_.last().label() << ": " <<
      new_w << " [" << algorithm_str_ << "]");
  if (have_expected_Lw_)
    Dout(dc::finish, " [expected_Lw: " << expected_Lw_ << "]");
  else
    Dout(dc::finish, "");
}
#endif

bool Algorithm::update_energy(double Lw)
{
  Dout(dc::warning, "Algorithm::update_energy: need to know if we need to check the energy");
  if (false/*state_ == IterationState::check_energy*/)
  {
    // Update the current kinetic energy. If this is an overshoot, abort this horizontal direction.
    //
    // This state is set when locally the curve looks like a cubic with its critical points
    // on the side that we just came from: in this case we move away from all critical points
    // and going uphill, losing energy. Therefore we need to check if we didn't go too high for
    // the kinetic energy that we have, in which case the exploration of this direction is aborted.

    // If the horizontal direction is still unknown, then we should always go towards the extreme
    // of the local approximation: in fact we would have jumped there and would not
    // care about the amount of available energy.
    ASSERT(hdirection_ != HorizontalDirection::undecided);

    // Update the energy and check if we had enough energy to reach this height.
    if (!energy_.maybe_update(Lw))
    {
      Dout(dc::notice, "Too much energy used: need to abort this direction (" << hdirection_ << ").");
      return false;
    }
  }
  else
  {
    // Update the current kinetic energy. If this is an overshoot, just increase the energy to 0.
    energy_.update(Lw);
  }
  return true;
}

bool Algorithm::handle_abort_hdirection(double& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_abort_hdirection(" << w << ")");

  // FIXME: implement.
  ASSERT(false);
  return false;
}

} // namespace gradient_descent
