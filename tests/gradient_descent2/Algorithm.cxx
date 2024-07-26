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
    event_server_.trigger(AlgorithmEventType{difference_event, w, expected_Lw_, Lw});
    have_expected_Lw_ = false;
  }
#endif

  // Create a new SampleNode for this sample.
  auto current = std::make_unique<SampleNode>(w, Lw, dLdw);
  Dout(dc::notice, "current = " << *current);

#ifdef CWDEBUG
  // Draw the new sample.
  event_server_.trigger(AlgorithmEventType{new_sample_event, *current});
#endif

  // Update kinetic energy. Returns false if too much energy was used.
  if (!update_energy(Lw))
    return handle_abort_hdirection(w);

  Dout(dc::notice, "Running " << state_);
  switch (state_)
  {
    case IterationState::initialization:
    {
      // We should only get here once.
      ASSERT(!last_);

      // Add the first sample.
      last_ = std::move(current);

#ifdef CWDEBUG
      {
        // Draw a line through this sample.
        math::CubicPolynomial cubic{Lw - dLdw * w, dLdw, 0.0, 0.0};
        event_server_.trigger(AlgorithmEventType{cubic_polynomial_event, cubic});
      }
#endif

      // Set w to the next value to probe. This uses last_.
      handle_single_sample(w);

      ASSERT(state_ == IterationState::next_sample);
      break;
    }

    case IterationState::next_sample:
    {
      last_ = std::move(current);
      w -= 1.0;
      expected_Lw_ = 1800.0;
      have_expected_Lw_ = true;
      Debug(set_algorithm_str(w, "minus one!"));
      break;
    }

    default:
      // Implement this state_.
      ASSERT(false);
  }

  return state_ != IterationState::success;
}

void Algorithm::handle_single_sample(double& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_single_sample(" << w << ")");

  double step;
#ifdef CWDEBUG
  char const* algorithm_str;
#endif

  if (small_step_ == 0.0)       // Not defined yet?
  {
    step = learning_rate_ * -last_->dLdw();

    // Did we drop on a (local) extreme as a starting point?!
    if (Scale::almost_zero(w, step))
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
    step = ((next_extreme_type_ == ExtremeType::minimum) == (last_->dLdw() > 0.0) ? -small_step_ : small_step_);
#ifdef CWDEBUG
    algorithm_str = "small step";
#endif
  }

#if 0 // FIXME: add scale back
  // This step could still be too small.
  if (last_->scale().negligible(step))
  {
    // This wouldn't be good because then the new sample will replace
    // the current one and we'd still have just one sample.
    step = last_->scale().make_significant(step);
#ifdef CWDEBUG
    algorithm_str = "avoiding replacement";
#endif
  }
#endif
  w += step;
  have_expected_Lw_ = false;
  Debug(set_algorithm_str(w, algorithm_str));
  state_ = IterationState::next_sample;
}

#ifdef CWDEBUG
void Algorithm::set_algorithm_str(double new_w, char const* algorithm_str)
{
  algorithm_str_ = algorithm_str;
  Dout(dc::notice|continued_cf, std::setprecision(12) << last_->w() << " --> " << last_->label() << ": " <<
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
