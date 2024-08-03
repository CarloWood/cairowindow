#include "sys.h"
#include "AnalyzedCubic.h"
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

      w = left_node->find_extreme(*right_node, next_extreme_type_);
      // w should be set if next_extreme_type_ was set.
      ASSERT(next_extreme_type_ == ExtremeType::unknown || w != SampleNode::uninitialized_magic);
      cubic_used_ = left_node;

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

        left_node->set_scale(scale_cp_type, critical_point_w, left_node->w(), right_node->w());
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
#ifdef CWDEBUG
      math::CubicPolynomial old_cubic = cubic_used_ != chain_.end() ? cubic_used_->cubic() : math::CubicPolynomial{};
#endif
      auto new_node = chain_.insert(std::move(current));
      auto right_node = std::next(new_node);

      AnalyzedCubic right_acubic;       // If right_node != chain_.end() is true then this corresponds to the cubic of std::next(new_node).
      bool right_has_extreme = false;
      double right_extreme_w;
      double right_dist;
      if (right_node != chain_.end())
      {
        new_node->initialize_cubic(*right_node COMMA_CWDEBUG_ONLY(event_server_, true));        // true: called on new_node.
        right_acubic.initialize(new_node->cubic(), next_extreme_type_);
        if (right_acubic.has_extremes())
        {
          right_extreme_w = right_acubic.get_extreme();
          right_has_extreme =
            (!left_of_ || right_extreme_w < left_of_->w()) &&
            (!right_of_ || right_of_->w() < right_extreme_w);
          right_dist = std::abs(right_extreme_w - new_node->w()) + std::abs(right_extreme_w - right_node->w());
        }
      }
      AnalyzedCubic left_acubic;        // If new_node != chain_.begin() is true then this corresponds to the cubic of std::prev(new_node).
      bool left_has_extreme = false;
      double left_extreme_w;
      double left_dist;
      if (new_node != chain_.begin())
      {
        auto left_node = std::prev(new_node);
        left_node->initialize_cubic(*new_node COMMA_CWDEBUG_ONLY(event_server_, false));        // false: new_node is passed.
        left_acubic.initialize(left_node->cubic(), next_extreme_type_);
        if (left_acubic.has_extremes())
        {
          left_extreme_w = left_acubic.get_extreme();
          left_has_extreme =
            (!left_of_ || left_extreme_w < left_of_->w()) &&
            (!right_of_ || right_of_->w() < left_extreme_w);
          left_dist = std::abs(left_extreme_w - left_node->w()) + std::abs(left_extreme_w - new_node->w());
        }
      }

      cubic_used_ = chain_.end();
      AnalyzedCubic const* used_acubic = nullptr;
      if (right_has_extreme)
      {
        // If both cubics have the required extreme in the required region, then
        // use the cubic whose sample points lay closest to their extreme.
        if (!left_has_extreme || right_dist < left_dist)
        {
          Dout(dc::notice|continued_cf, "Choosing " << next_extreme_type_ << " of right cubic because ");
#ifdef CWDEBUG
          if (left_has_extreme)
            Dout(dc::finish, "that extreme (" << right_extreme_w << ") is closer to the center of [" <<
                new_node->label() << "]<--->[" << right_node->label() << "], than the left cubic extreme (" <<
                left_extreme_w << ") is from the center of [" << std::prev(new_node)->label() << "]<--->[" << new_node->label() <<
                "]; (" << right_dist << " < " << left_dist << ").");
          else if (new_node != chain_.begin())
            Dout(dc::finish, "the left one does not have that extreme between " <<
                (right_of_ ? right_of_->w() : -std::numeric_limits<double>::infinity()) << " and " <<
                (left_of_ ? left_of_->w() : std::numeric_limits<double>::infinity()));
          else
            Dout(dc::finish, "the left one does not exists.");
#endif
          w = right_extreme_w;
          cubic_used_ = new_node;
          used_acubic = &right_acubic;
        }
        else
        {
          Dout(dc::notice|continued_cf, "Choosing " << next_extreme_type_ << " of left cubic because ");
#ifdef CWDEBUG
          if (left_has_extreme)
            Dout(dc::finish, "that extreme (" << left_extreme_w << ") is closer to the center of [" <<
                std::prev(new_node)->label() << "]<--->[" << new_node->label() << "], than the right cubic extreme (" <<
                right_extreme_w << ") is from the center of " << new_node->label() << "]<--->[" << right_node->label() <<
                "; (" << left_dist << " < " << right_dist << ").");
          else
            Dout(dc::finish, "the right one does not have that extreme between " <<
                (right_of_ ? right_of_->w() : -std::numeric_limits<double>::infinity()) << " and " <<
                (left_of_ ? left_of_->w() : std::numeric_limits<double>::infinity()));
#endif
          w = left_extreme_w;
          cubic_used_ = std::prev(new_node);
          used_acubic = &left_acubic;
        }
      }
      else
      {
        if (left_has_extreme)
        {
#ifdef CWDEBUG
          Dout(dc::notice|continued_cf, "Choosing " << next_extreme_type_ << " of left cubic because the right one does not ");
          if (right_node != chain_.end())
            Dout(dc::finish, "have that extreme between " <<
              (right_of_ ? right_of_->w() : -std::numeric_limits<double>::infinity()) << " and " <<
              (left_of_ ? left_of_->w() : std::numeric_limits<double>::infinity()));
          else
            Dout(dc::finish, "exists.");
#endif
          w = left_extreme_w;
          cubic_used_ = std::prev(new_node);
          used_acubic = &left_acubic;
        }
        else
        {
          ASSERT(false);
          break;
        }
      }
      if (w < new_node->w())
        left_of_ = &*new_node;
      else
        right_of_ = &*new_node;

#ifdef CWDEBUG
      event_server_.trigger(AlgorithmEventType{left_of_right_of_event, left_of_, right_of_,
          next_extreme_type_, w, cubic_used_->cubic()(w)});
#endif

      // Update scale.
      if (cubic_used_ != chain_.end())
      {
        CriticalPointType scale_cp_type = next_extreme_type_ == ExtremeType::minimum ? CriticalPointType::minimum
                                                                                     : CriticalPointType::maximum;
        math::CubicPolynomial const& cubic = cubic_used_->cubic();

        // Find the left/right edge sample as far away from the samples used to fit the current cubic, that still matches that cubic.
        //
        // Start with setting left_edge to cubic_used_ and right_edge to the next sample:
        // respectively left and right sample that were used to fit the cubic.
        //
        //        \      /
        //         \    🞄
        //          `-🞄´^
        //            ^ right
        //         left
        SampleNode::const_iterator left_edge = cubic_used_;
        SampleNode::const_iterator right_edge = std::next(cubic_used_);

#ifdef CWDEBUG
        // Because these points were used to generate the cubic, they should perfectly fit it.
        if (!Scale::matches(left_edge->w(), left_edge->Lw(), cubic, *used_acubic))
        {
          Dout(dc::warning, "left_edge, " << *left_edge << " is not matching " << cubic << " / " << *used_acubic);
        }
        if (!Scale::matches(right_edge->w(), right_edge->Lw(), cubic, *used_acubic))
        {
          Dout(dc::warning, "right_edge, " << *right_edge << " is not matching " << cubic << " / " << *used_acubic);
        }
#endif
        // Now advance left_edge to the left as far as possible, but not
        // so far that it points to a sample that does not match.
        while (left_edge != chain_.begin())
        {
          --left_edge;
          if (!Scale::matches(left_edge->w(), left_edge->Lw(), cubic, *used_acubic))
          {
            ++left_edge;
            break;
          }
        }
        // Same with right_edge, but as far as possible to the right.
        while (++right_edge != chain_.end() && Scale::matches(right_edge->w(), right_edge->Lw(), cubic, *used_acubic))
          ;
        // In this case we did overshoot by one, so move right_edge one back.
        auto next_right_edge = right_edge--;

        // Thus, left_edge and right_edge as the last "good" (matching) samples.
        double left_edge_w_good = left_edge->w();
        double right_edge_w_good = right_edge->w();

        // Is there a bad left sample?
        if (left_edge != chain_.begin())
        {
          auto left_bad_sample = std::prev(left_edge);
          double left_edge_w = left_bad_sample->w();
          double left_edge_Lw;
          do
          {
            // Half the distance to left_edge_w until it matches again.
            left_edge_w = 0.5 * (left_edge_w_good + left_edge_w);
            left_edge_Lw = left_bad_sample->cubic()(left_edge_w);
          }
          while (!Scale::matches(left_edge_w, left_edge_Lw, cubic, *used_acubic));
          // Use this as the left edge.
          left_edge_w_good = left_edge_w;
        }
        // Is there a bad right sample?
        if (next_right_edge != chain_.end())
        {
          double right_edge_w = next_right_edge->w();   // next_right_edge is the bad sample.
          double right_edge_Lw;
          do
          {
            right_edge_w = 0.5 * (right_edge_w_good + right_edge_w);
            right_edge_Lw = right_edge->cubic()(right_edge_w);
          }
          while (!Scale::matches(right_edge_w, right_edge_Lw, cubic, *used_acubic));
          // Use this as the right edge.
          right_edge_w_good = right_edge_w;
        }

        cubic_used_->set_scale(scale_cp_type, w, left_edge_w_good, right_edge_w_good);
#ifdef CWDEBUG
        event_server_.trigger(AlgorithmEventType{scale_draw_event, ScaleUpdate::towards_cp, *cubic_used_, old_cubic});
#endif
      }

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
