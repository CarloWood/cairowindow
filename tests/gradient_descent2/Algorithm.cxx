#include "sys.h"
#include "AnalyzedCubic.h"
#include "SampleNode.h"
#include "Algorithm.h"
#include "IterationState.h"
#include <Eigen/Dense>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/at_scope_end.h"
#include "cwds/Restart.h"
#endif

namespace gradient_descent {

bool Algorithm::operator()(double& w, double Lw, double dLdw)
{
  DoutEntering(dc::notice, "Algorithm::operator()(" << w << ", " << Lw << ", " << dLdw << ")");
  RESTART

  // Do not re-enter this function after it returned false!
  ASSERT(state_ != IterationState::success);

#ifdef CWDEBUG
  // Erase all previous curves (if they exist).
  event_server_.trigger(AlgorithmEventType{reset_event});
  if (have_expected_Lw_)
  {
    // Draw an arrow from where expected to end up, based on the approximation, and the actual value of L(w).
    event_server_.trigger(AlgorithmEventType{difference_event, w, expected_Lw_, Lw});
    have_expected_Lw_ = false;
  }

  // Print the chain_ to debug output and do a sanity check just before leaving this function.
  auto&& dump_chain = at_scope_end([this]{ chain_.dump(this); chain_.sanity_check(this); });
#endif

  // Create a new Sample for this sample.
  Sample current{w, Lw, dLdw COMMA_CWDEBUG_ONLY(label_context_)};
  Dout(dc::notice, "current = " << current);

#ifdef CWDEBUG
  // Draw the new sample.
  event_server_.trigger(AlgorithmEventType{new_sample_event, current});
#endif

  if (state_ == IterationState::initialization)
  {
    // We should only get here the first time.
    ASSERT(chain_.empty());

    DEBUG_ONLY(bool success =) update_energy(Lw);
    // Can never want to abort on energy when still initializing.
    ASSERT(success);

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
    chain_.find_larger(w);

    state_ = IterationState::first_cubic;
    return true;
  }

  // Insert the new sample.
  auto new_node = chain_.insert(std::move(current));

  // Initialize the cubic(s).
  SampleNode::const_iterator right_node = std::next(new_node);  // The node right of the new node.
  SampleNode::const_iterator left_node;                         // The node left of the new node (only valid if new_node != chain_.begin()).
  if (new_node != chain_.begin())
  {
    left_node = std::prev(new_node);
    left_node->initialize_cubic(new_node COMMA_CWDEBUG_ONLY(event_server_, false));        // false: new_node is passed.
  }
  if (right_node != chain_.end())
    new_node->initialize_cubic(right_node COMMA_CWDEBUG_ONLY(event_server_, true));        // true: called on new_node.

  for (;;)
  {
    // Update kinetic energy. Returns false if too much energy was used.
    if (!update_energy(Lw))
    {
      if (!handle_abort_hdirection(w))
        return false;
    }
    else
    {
      Dout(dc::notice, "⁂  Running \"" << state_ << '"');
      switch (state_)
      {
        case IterationState::first_cubic:
        {
          auto cubic_used_ = right_node == chain_.end() ? left_node : new_node;
          auto next = std::next(cubic_used_);

          // Jump to extreme of the first cubic.
          w = cubic_used_->find_extreme(*next, next_extreme_type_);

          // w should be set if next_extreme_type_ was set.
          ASSERT(next_extreme_type_ == ExtremeType::unknown || w != SampleNode::uninitialized_magic);

          // Initialize the scale.
          {
            CriticalPointType scale_cp_type =
              next_extreme_type_ == ExtremeType::unknown ? CriticalPointType::inflection_point
                                                         : next_extreme_type_ == ExtremeType::minimum ? CriticalPointType::minimum
                                                                                                      : CriticalPointType::maximum;
            double critical_point_w =
              next_extreme_type_ == ExtremeType::unknown ? cubic_used_->cubic().inflection_point()
                                                         : w;

            cubic_used_->set_scale(scale_cp_type, critical_point_w, cubic_used_->w(), next->w());
#ifdef CWDEBUG
            event_server_.trigger(AlgorithmEventType{scale_draw_event, ScaleUpdate::initialized, *cubic_used_, {}});
#endif
          }

          // If cubic_used_ didn't have a usable extreme, then next_extreme_type_ was set to unknown.
          if (next_extreme_type_ == ExtremeType::unknown)
          {
            // We're going to look for a minimum.
            if (cubic_used_->is_rising())
            {
              // Keep going to the left.
              w = cubic_used_->w() - cubic_used_->scale().value();
              left_of_ = cubic_used_;
            }
            else // cubic_used_->is_falling() or the cubic is flat.
            {
              // Keep going to the right.
              w = next->w() + cubic_used_->scale().value();
              if (cubic_used_->is_falling())
                right_of_ = next;
            }
          }
          else
          {
            ASSERT(w != SampleNode::uninitialized_magic);
            double const negligible_offset = negligible_scale_fraction * cubic_used_->scale().value();
            if (cubic_used_->w() - w > negligible_offset)
              left_of_ = cubic_used_;
            else if (next->w() - w > negligible_offset)
            {
              left_of_ = next;
              if (w - cubic_used_->w() > negligible_offset)
                right_of_ = cubic_used_;
            }
            else if (w - next->w() > negligible_offset)
              right_of_ = next;
          }
#ifdef CWDEBUG
          event_server_.trigger(AlgorithmEventType{left_of_right_of_event,
              left_of_ == chain_.end() ? nullptr : &*left_of_,
              right_of_ == chain_.end() ? nullptr : &*right_of_,
              next_extreme_type_, w, cubic_used_->cubic()(w)});
#endif

          // If next_extreme_type_ was not set then this cubic has no extremes. In that case go for a minimum.
          if (next_extreme_type_ == ExtremeType::unknown)
            next_extreme_type_ = ExtremeType::minimum;

          state_ = IterationState::find_extreme;

          // This means we already found a local extreme.
          if (std::abs(w - new_node->w()) <= 1e-3 * cubic_used_->scale().value())     // Too close to the last sample?
          {
            [[maybe_unused]] bool success = handle_local_extreme(w);
            // We only added two samples, so handle_local_extreme is just going to ask for a another sample.
            ASSERT(success);
          }

          break;
        }

        case IterationState::find_extreme:
        {
#ifdef CWDEBUG
          math::CubicPolynomial old_cubic = cubic_used_ != chain_.end() ? cubic_used_->cubic() : math::CubicPolynomial{};
#endif
          AnalyzedCubic right_acubic;       // If right_node != chain_.end() is true then this corresponds to the cubic of std::next(new_node).
          bool right_has_extreme = false;
          double right_extreme_w;
          double right_dist;
          if (right_node != chain_.end())
          {
            right_acubic.initialize(new_node->cubic(), next_extreme_type_);
            if (right_acubic.has_extremes())
            {
              //FIXME: initialize scale immediately after an initialize_cubic? And then use the real scale here.
              double const significant_offset = significant_scale_fraction * new_node->w_scale_estimate();
              right_extreme_w = right_acubic.get_extreme();
              right_has_extreme =
                (left_of_ == chain_.end() || right_extreme_w - left_of_->w() < significant_offset) &&
                (right_of_ == chain_.end() || right_of_->w() - right_extreme_w < significant_offset);
              right_dist = std::abs(right_extreme_w - new_node->w()) + std::abs(right_extreme_w - right_node->w());
            }
          }
          AnalyzedCubic left_acubic;        // If new_node != chain_.begin() is true then this corresponds to the cubic of std::prev(new_node).
          bool left_has_extreme = false;
          double left_extreme_w;
          double left_dist;
          if (new_node != chain_.begin())
          {
            left_acubic.initialize(left_node->cubic(), next_extreme_type_);
            if (left_acubic.has_extremes())
            {
              //FIXME: initialize scale immediately after an initialize_cubic? And then use the real scale here.
              double const significant_offset = significant_scale_fraction * left_node->w_scale_estimate();
              left_extreme_w = left_acubic.get_extreme();
              left_has_extreme =
                (left_of_ == chain_.end() || left_extreme_w - left_of_->w() < significant_offset) &&
                (right_of_ == chain_.end() || right_of_->w() - left_extreme_w < significant_offset);
              left_dist = std::abs(left_extreme_w - left_node->w()) + std::abs(left_extreme_w - new_node->w());
            }
          }

          cubic_used_ = chain_.end();
          AnalyzedCubic* used_acubic = nullptr;
          HorizontalDirection keep_going = HorizontalDirection::undecided;      // Default: undecided means nothing needs to be done.
          if (right_has_extreme)
          {
            // If both cubics have the required extreme in the required region, then
            // use the cubic whose samples lay closest to their extreme.
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
                    (right_of_ != chain_.end() ? right_of_->w() : -std::numeric_limits<double>::infinity()) << " and " <<
                    (left_of_ != chain_.end() ? left_of_->w() : std::numeric_limits<double>::infinity()));
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
                    right_extreme_w << ") is from the center of [" << new_node->label() << "]<--->[" << right_node->label() <<
                    "]; (" << left_dist << " < " << right_dist << ").");
              else
                Dout(dc::finish, "the right one does not have that extreme between " <<
                    (right_of_ != chain_.end() ? right_of_->w() : -std::numeric_limits<double>::infinity()) << " and " <<
                    (left_of_ != chain_.end() ? left_of_->w() : std::numeric_limits<double>::infinity()));
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
                  (right_of_ != chain_.end() ? right_of_->w() : -std::numeric_limits<double>::infinity()) << " and " <<
                  (left_of_ != chain_.end() ? left_of_->w() : std::numeric_limits<double>::infinity()));
              else
                Dout(dc::finish, "exists.");
#endif
              w = left_extreme_w;
              cubic_used_ = std::prev(new_node);
              used_acubic = &left_acubic;
            }
            else
            {
#ifdef CWDEBUG
              if (new_node == chain_.begin())
              {
                ASSERT(right_node != chain_.end());
                Dout(dc::notice|continued_cf, "There is no left cubic, and the right cubic has no extreme between ");
              }
              else if (right_node == chain_.end())
                Dout(dc::notice|continued_cf, "There is no right cubic, and the left cubic has no extreme between ");
              else
                Dout(dc::notice|continued_cf, "Neither the left nor the right cubic has an extreme between ");
              Dout(dc::finish, (right_of_ != chain_.end() ? right_of_->w() : -std::numeric_limits<double>::infinity()) <<
                  " and " << (left_of_ != chain_.end() ? left_of_->w() : std::numeric_limits<double>::infinity()));
#endif
              auto current = new_node;
              // If we are on the wrong side of 'foo_of_', jump to foo_of_.
              //
              //                           |==>           <==|
              //           ------x------right_of-----x-----left_of------x-------
              //                 |--jump-->|       stay      |<---jump--|
              if (right_of_ != chain_.end() && current->w() < right_of_->w())
                current = right_of_;
              else if (left_of_ != chain_.end() && current->w() > left_of_->w())
                current = left_of_;

              // Determine the direction that we must go in. This must rely on the type_ of the cubics
              // rather than the sign of the derivative because the latter can't be trusted to be consistent
              // (for example, it could be -1e-14; which is close enough to zero that it might have ended up
              // less than zero dure to floating-point round-off errors).

              HorizontalDirection direction;
              if (left_of_ == chain_.end() && right_of_ == chain_.end())
              {
                // If no region boundaries are set, then the cubics don't have an extreme at all; aka,
                // we are in the middle of a monotonic rising or falling area: in that case just
                // move the current point to the highest/lowest point left/right of the current sample,
                // depending on whether we're looking for a minimum or maximum.
                bool at_right_edge = current->type() == CubicToNextSampleType::unknown;
                if (at_right_edge)
                  current = std::prev(current);
                ASSERT(current->type() != CubicToNextSampleType::unknown);
                direction =
                  (next_extreme_type_ == ExtremeType::minimum && current->is_falling()) ||
                  (next_extreme_type_ == ExtremeType::maximum && current->is_rising()) ?
                  HorizontalDirection::right : HorizontalDirection::left;
                if (at_right_edge)
                  current = std::next(current);
              }
              else
              {
                // If there is one region boundary, then we should start from that boundary and move
                // away from it until we find a cubic with the sought for local extreme, or reach the
                // edge on the other side - in which case we'll have to "keep going" again by adding
                // the scale (and return to probe the next sample).
                if (left_of_ == chain_.end())
                  current = right_of_;
                else if (right_of_ == chain_.end())
                  current = left_of_;

                // It should never happen that both left_of_ and right_of_ are set,
                // while we're not already on one of those boundaries.
                ASSERT(current == left_of_ || current == right_of_);

                // In this case we just want to move away from the boundary that we're on.
                direction = current == right_of_ ? HorizontalDirection::right : HorizontalDirection::left;
              }
              Dout(dc::notice, "direction = " << direction);

              // Move current in the given direction until we find a cubic with the sought for local extreme, or reach the end of the chain.
              if (direction == HorizontalDirection::right)
              {
                while (current->type() != CubicToNextSampleType::unknown && !current->has_extreme(next_extreme_type_))
                  current = std::next(current);
                if (current->type() == CubicToNextSampleType::unknown)
                {
                  right_of_ = current;
                  // Set cubic_used_ in order to set the scale on it: if we keep going again,
                  // we want to use the scale of this cubic.
                  cubic_used_ = std::prev(current);
                  right_acubic.initialize(cubic_used_->cubic(), next_extreme_type_);
                  used_acubic = &right_acubic;
                  // Keep going to the right.
                  w = right_of_->w();
                  keep_going = HorizontalDirection::right;      // cubic_used_->scale().value() needs to be added to w.
                }
                else
                {
                  // implement
                  ASSERT(false);
                }
              }
              else
              {
                while (current != chain_.begin() && !std::prev(current)->has_extreme(next_extreme_type_))
                  current = std::prev(current);
                if (current == chain_.begin())
                {
                  left_of_ = current;
                  // Set cubic_used_ in order to set the scale on it: if we keep going again,
                  // we want to use the scale of this cubic.
                  cubic_used_ = current;
                  left_acubic.initialize(cubic_used_->cubic(), next_extreme_type_);
                  used_acubic = &left_acubic;
                  // Keep going to the left.
                  w = left_of_->w();
                  keep_going = HorizontalDirection::left;       // cubic_used_->scale().value() needs to be subtracted from w.
                }
                else
                {
                  // implement
                  ASSERT(false);
                }
              }
            }
          }

          // cubic_used_ should have been set to the cubic that we need to update the scale of.
          ASSERT(cubic_used_ != chain_.end());

          // Update scale.
          CriticalPointType const scale_cp_type =
            used_acubic->has_extremes() ? (next_extreme_type_ == ExtremeType::minimum ? CriticalPointType::minimum
                                                                                      : CriticalPointType::maximum)
                                        : CriticalPointType::inflection_point;

          double const critical_point_w =
            used_acubic->has_extremes() ? used_acubic->get_extreme()
                                        : used_acubic->inflection_point();

          math::CubicPolynomial const& cubic = cubic_used_->cubic();

          // Find the left/right edge sample as far away from the samples used to fit the current cubic, that still matches that cubic.
          //
          // Start with setting left_edge to cubic_used_ and right_edge to the next sample,
          // respectively the left and right sample that were used to fit the cubic.
          //
          //        \      /
          //         \    🞄
          //          `🞄-´^
          //           ^  right
          //        left
          SampleNode::const_iterator left_edge = cubic_used_;
          SampleNode::const_iterator right_edge = std::next(cubic_used_);

          // Must be called before calling matches.
          used_acubic->initialize_matches(*left_edge, *right_edge);

#ifdef CWDEBUG
          // Because these points were used to generate the cubic, they should perfectly fit it.
          if (!used_acubic->matches(left_edge->w(), left_edge->Lw(), cubic))
          {
            Dout(dc::warning, "left_edge, " << *left_edge << " is not matching " << cubic << " / " << *used_acubic);
          }
          if (!used_acubic->matches(right_edge->w(), right_edge->Lw(), cubic))
          {
            Dout(dc::warning, "right_edge, " << *right_edge << " is not matching " << cubic << " / " << *used_acubic);
          }
#endif
          // Now advance left_edge to the left as far as possible, but not
          // so far that it points to a sample that does not match.
          while (left_edge != chain_.begin())
          {
            --left_edge;
            if (!used_acubic->matches(left_edge->w(), left_edge->Lw(), cubic))
            {
              ++left_edge;
              break;
            }
          }
          // Same with right_edge, but as far as possible to the right.
          while (++right_edge != chain_.end() && used_acubic->matches(right_edge->w(), right_edge->Lw(), cubic))
            ;
          // In this case we did overshoot by one, so move right_edge one back.
          auto next_right_edge = right_edge--;

          // Thus, left_edge and right_edge are the last "good" (matching) samples.
          double left_edge_w_good = left_edge->w();
          double right_edge_w_good = right_edge->w();

          // Is there a (bad) sample on the left?
          if (left_edge != chain_.begin())
          {
            auto left_bad_sample = std::prev(left_edge);
            double left_edge_w = left_bad_sample->w();
            double left_edge_Lw;
            int guard = 32;
            do
            {
              // Half the distance to left_edge_w until it matches again.
              left_edge_w = 0.5 * (left_edge_w_good + left_edge_w);
              // Eventually left_edge_w becomes equal to left_edge_w_good; however, due to
              // floating-point round-off errors even matches(left_edge_w_good, cubic)
              // might be false; therefore stop regardless after 32 loops.
              if (--guard == 0) break;
              left_edge_Lw = left_bad_sample->cubic()(left_edge_w);
            }
            while (!used_acubic->matches(left_edge_w, left_edge_Lw, cubic));
            // Use this as the left edge.
            left_edge_w_good = left_edge_w;
          }
          // Is there a (bad) sample on the right?
          if (next_right_edge != chain_.end())
          {
            double right_edge_w = next_right_edge->w();   // next_right_edge is the bad sample.
            double right_edge_Lw;
            int guard = 32;
            do
            {
              right_edge_w = 0.5 * (right_edge_w_good + right_edge_w);
              if (--guard == 0) break;
              right_edge_Lw = right_edge->cubic()(right_edge_w);
            }
            while (!used_acubic->matches(right_edge_w, right_edge_Lw, cubic));
            // Use this as the right edge.
            right_edge_w_good = right_edge_w;
          }

          cubic_used_->set_scale(scale_cp_type, critical_point_w, left_edge_w_good, right_edge_w_good);
#ifdef CWDEBUG
          event_server_.trigger(AlgorithmEventType{scale_draw_event, ScaleUpdate::towards_cp, *cubic_used_, old_cubic});
#endif

          // Add scale if must "keep going" (there was no extreme to jump to).
          w += static_cast<int>(keep_going) * cubic_used_->scale().value();

          // Adjust left_of_/right_of_.
          double const negligible_offset = negligible_scale_fraction * cubic_used_->scale().value();
          if (new_node->w() - w > negligible_offset)
            left_of_ = new_node;
          else if (w - new_node->w() > negligible_offset)
            right_of_ = new_node;

#ifdef CWDEBUG
          event_server_.trigger(AlgorithmEventType{left_of_right_of_event,
              left_of_ == chain_.end() ? nullptr : &*left_of_,
              right_of_ == chain_.end() ? nullptr : &*right_of_,
              next_extreme_type_, w, cubic_used_->cubic()(w)});
#endif

          // If the last step is significantly smaller than the scale, then we found a local extreme.
          double step = std::abs(chain_.last()->w() - w);
          Dout(dc::notice, "step = " << step);
          if (step < (next_extreme_type_ == ExtremeType::minimum ? 0.01 : 0.05) * cubic_used_->scale().value() && !handle_local_extreme(w) && !handle_abort_hdirection(w))
            return false;

          break;
        }

        case IterationState::extra_sample:
        {
          if (!handle_local_extreme(w) && !handle_abort_hdirection(w))
            return false;

          break;
        }

        case IterationState::finish:
        {
          // Store the last sample and set the state to success.
          set_global_minimum(new_node);
          return false;
        }

        default:
          // Implement this state_.
          ASSERT(false);
      }
    }

    // Already find the point where to insert the next sample (with value w).
    // We do this by looking for the first SampleNode that has a w value that
    // is greater than that of the new sample and then inserting current in
    // front of that. This function only finds that larger sample.
    chain_.find_larger(w);

    // Check if this is a duplicate (closer than 0.00001 * scale) to an existing sample.
    auto ibp = chain_.duplicate(cubic_used_->scale().value());
    if (AI_LIKELY(!ibp.second) || (state_ == IterationState::finish && ibp.first->is_fake()))
    {
      // Ask for a new sample at `w` if not already finished.
      Dout(dc::notice, "Next probe: " << std::setprecision(std::numeric_limits<long double>::max_digits10) << w);
      break;
    }

    Dout(dc::notice, "New probe (" << w << ") is too close to an existing sample (" << *ibp.first << "), reusing that.");
    new_node = ibp.first;
    chain_.reuse(new_node);
  }

  return true;
}

void Algorithm::handle_single_sample(double& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_single_sample(" << w << ")");

  Sample const& sample = *chain_.last();
  double step;
#ifdef CWDEBUG
  char const* algorithm_str;
#endif

  if (small_step_ == 0.0)       // Not defined yet?
  {
    step = -learning_rate_ * sample.dLdw();

    // Did we drop on a (local) extreme as a starting point?!
    if (std::abs(step) < epsilon || std::abs(step) < negligible_scale_fraction * std::abs(w))
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
  Dout(dc::notice|continued_cf, std::setprecision(12) << chain_.last()->w() << " --> " << chain_.last()->label() << ": " <<
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

  // Has a minimum been found at all?
  if (best_minimum_ == chain_.end())
  {
    Dout(dc::notice, "Aborting going " << hdirection_ << " and terminating search because no minimum has has been found at all!");
    return false;
  }

  // Jump back to the best minimum and continue in the opposite hdirection.
  Dout(dc::notice, "Aborting exploring " << hdirection_ << " of the minimum at " << static_cast<Sample const&>(*best_minimum_) << ".");

  // Was this minimum already explored in both directions?
  if (best_minimum_->done())
    return false;

  ASSERT(hdirection_ != HorizontalDirection::undecided);

  // Change hdirection.
  hdirection_ = opposite(hdirection_);
  Dout(dc::notice, "Changed horizontal direction to " << hdirection_);
  last_extreme_ = best_minimum_;
  cubic_used_ = chain_.end();
  left_of_ = chain_.end();
  right_of_ = chain_.end();
  // Next we'll be looking for a maximum.
  next_extreme_type_ = ExtremeType::maximum;
  //small_step_ = 

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{left_of_right_of_event,
      left_of_ == chain_.end() ? nullptr : &*left_of_,
      right_of_ == chain_.end() ? nullptr : &*right_of_,
      next_extreme_type_, w, cubic_used_->cubic()(w)});
  event_server_.trigger(AlgorithmEventType{hdirection_known_event, *last_extreme_, hdirection_});
#endif

  return true;
}

bool Algorithm::handle_local_extreme(double& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_local_extreme(" << w << ")" << (state_ == IterationState::extra_sample ? " [Extra Sample]" : ""));
  RESTART

  double const scale_value = cubic_used_->scale().value();

  auto extra_sample = state_ != IterationState::extra_sample ? chain_.end() : chain_.last();

  if (state_ != IterationState::extra_sample)
  {
    // This block sets small_step_, last_extreme_ and best_minimum_.

    // Update small_step_ (value() returns a positive value).
    small_step_ = scale_value;
    Dout(dc::notice, "small_step_ set to " << small_step_);

    // Set chain_.larger_ to point to the first SampleNode that is larger than w, if any.
    chain_.find_larger(w);

    // Find out if a new (fake) SampleNode should inserted.
    SampleNode::const_iterator right = chain_.larger();
    SampleNode::const_iterator left = right == chain_.begin() ? chain_.end() : std::prev(right);

    double const negligible_offset = negligible_scale_fraction * scale_value;

    if (right != chain_.end() && std::abs(right->w() - w) <= negligible_offset)
    {
      Dout(dc::notice, "Using existing sample on the right at " << right->w());
      last_extreme_ = right;
    }
    else if (left != chain_.end() && std::abs(left->w() - w) <= negligible_offset)
    {
      Dout(dc::notice, "Using existing sample on the left at " << left->w());
      last_extreme_ = left;
    }
    else
    {
      // Create a new Sample for this sample.
      Sample fake_sample{w, cubic_used_->cubic()(w), 0.0 COMMA_CWDEBUG_ONLY(label_context_)};
      Dout(dc::notice, "Using fake sample = " << fake_sample);
      last_extreme_ = chain_.insert(std::move(fake_sample));
      last_extreme_->set_fake(true);
      auto next = std::next(last_extreme_);
      if (next != chain_.end())
        last_extreme_->initialize_cubic(next, event_server_, true);
      if (last_extreme_ != chain_.begin())
        std::prev(last_extreme_)->initialize_cubic(last_extreme_ COMMA_CWDEBUG_ONLY(event_server_, false));
    }

    // Mark this sample as local extreme.
    last_extreme_->set_local_extreme(next_extreme_type_);

#ifdef CWDEBUG
    // Draw the local extreme sample.
    event_server_.trigger(AlgorithmEventType{new_local_extreme_event, *last_extreme_, std::to_string(last_extreme_->label())});
    event_server_.trigger(AlgorithmEventType{scale_erase_event});

    // Do some sanity checks.
    if (next_extreme_type_ == ExtremeType::minimum)
    {
      // This is a local minimum.
      ASSERT(last_extreme_->type() == CubicToNextSampleType::unknown || last_extreme_->has_right_min());
      ASSERT(last_extreme_ == chain_.begin() || std::prev(last_extreme_)->has_left_min());

      // We always want a scale to be defined for a minimum.
      ASSERT(std::next(last_extreme_) == chain_.end() || last_extreme_->scale().type() != CriticalPointType::none);
    }
    else
    {
      // This is a local maximum.
      ASSERT(next_extreme_type_ == ExtremeType::maximum);
      ASSERT(last_extreme_->type() == CubicToNextSampleType::unknown || last_extreme_->has_right_max());
      ASSERT(last_extreme_ == chain_.begin() || std::prev(last_extreme_)->has_left_max());
    }
#endif

    // Keep track of the best minimum so far; or abort if this minimum isn't better than one found before.
    if (next_extreme_type_ == ExtremeType::minimum)
    {
      Dout(dc::notice, "new extreme = " << *last_extreme_);
      // We were looking for a minimum.
      ASSERT(last_extreme_->get_extreme_type() == ExtremeType::minimum);
      if (best_minimum_ == chain_.end() || best_minimum_->Lw() > last_extreme_->Lw())
      {
        best_minimum_ = last_extreme_;
        Dout(dc::notice, "best_minimum_ set to " << *best_minimum_ << " and cubic approximation: " << *cubic_used_);
      }
      if (last_extreme_ != best_minimum_)
      {
        // The new minimum isn't better than what we found already. Stop going into this direction.
        Dout(dc::notice, "The new minimum (at " << static_cast<Sample const&>(*last_extreme_) << ") isn't better than what we found already. "
            "Stop going into the direction " << hdirection_ << ".");
        return false;
      }
    }
  }
#ifdef CWDEBUG
  else
    // The w value is still at that of the (requested) extra sample.
    ASSERT(extra_sample->w() == w);
#endif

  // Find all samples around last_extreme_ that are within scale.
  SampleNode::const_iterator left_edge = last_extreme_;
  SampleNode::const_iterator right_edge = last_extreme_;
  double const extreme_w = last_extreme_->w();

  static int constexpr max_samples_per_side = 3;
  std::array<SampleNode::const_iterator, 2 * max_samples_per_side + 1> usable_samples;
  int max_size = 2 * max_samples_per_side;

  // Always add the found local extreme, unless it is a fake point.
  int one_passed_extreme_index = max_samples_per_side;          // Or max_samples_per_side if the last_extreme_ isn't stored because it is fake.
  if (!last_extreme_->is_fake())
  {
    usable_samples[one_passed_extreme_index++] = last_extreme_;
    ++max_size;
  }

  // Now advance left_edge to the left as far as possible, but not
  // so far that it points to a sample that does not match.
  double const left_scale = cubic_used_->scale().left_edge_w();       // The exact w() value of the left-most Sample used for the scale.
  int left_i = max_samples_per_side;
  while (left_edge != chain_.begin())
  {
    --left_edge;
    // Skip if the sample is too far to the left and is not the extra sample that we just requested.
    if (left_edge->w() < left_scale && left_edge != extra_sample)
    {
      ++left_edge;
      break;
    }
    // Calculate the distance to the previous sample, if any.
    double dist = left_i == max_samples_per_side ? scale_value /*large enough*/ :  std::abs(left_edge->w() - usable_samples[left_i]->w());
    // If the current sample is fake or too close to the previous sample, then skip this sample.
    if (!left_edge->is_fake() && dist >= 0.001 * scale_value)
    {
      if (left_i == 0)
        break;          // Overflow
      usable_samples[--left_i] = left_edge;
    }
  }
  // Same with right_edge, but as far as possible to the right.
  double const right_scale = cubic_used_->scale().right_edge_w();       // The exact w() value of the right-most Sample used for the scale.
  int right_i = one_passed_extreme_index;
  while (++right_edge != chain_.end() && (right_edge->w() <= right_scale || right_edge == extra_sample))
  {
    // Calculate the distance to the previous sample, if any.
    double dist = (right_i == one_passed_extreme_index && left_i == max_samples_per_side) ?
      scale_value /*large enough*/ : std::abs(right_edge->w() - usable_samples[right_i - 1]->w());
    // If the current sample is fake or too close to the previous sample, then skip this sample.
    if (!right_edge->is_fake() && dist >= 0.001 * scale_value)
    {
      if (right_i == max_size)
        break;          // Overflow
      usable_samples[right_i++] = right_edge;
    }
  }
  // In this case we did overshoot by one, so right_edge points passed the last useable sample.

  int number_of_usable_samples = right_i - left_i;
#ifdef CWDEBUG
  Dout(dc::notice, "Number of samples within scale range: " << number_of_usable_samples);
  for (int i = left_i; i < right_i; ++i)
    Dout(dc::notice, *usable_samples[i]);
#endif

  // There must be at least one usable (real) sample: one of the samples used to fit cubic_used_ may lay further away than 1.1 times the scale, but not both.
  ASSERT(number_of_usable_samples > 0);

  // Set local_extreme_index such that usable_samples[local_extreme_index] is the sample that is the closest to local_extreme_.
  // If last_extreme_->is_fake() is false then last_extreme_ itself is stored at index max_samples_per_side.
  // If left_i == max_samples_per_side then no 'left' sample was stored, and the closest right sample was stored at index max_samples_per_side.
  // If right_i == max_samples_per_side then no 'right' sample was stored, and the closest left sample was stored at index max_samples_per_side - 1.
  int const local_extreme_index = (last_extreme_->is_fake() && left_i < max_samples_per_side &&
      (right_i == max_samples_per_side || std::abs(usable_samples[max_samples_per_side - 1]->w() - extreme_w) < std::abs(usable_samples[max_samples_per_side]->w() - extreme_w))) ?
    max_samples_per_side - 1 : max_samples_per_side;

  // After asking for one extra sample we should now have two or three samples.
  ASSERT(state_ != IterationState::extra_sample || number_of_usable_samples >= 2);

  // If we do not already have at least three relevant samples then get another one.
  if (number_of_usable_samples < 3)
  {
    ASSERT(number_of_usable_samples == 2);

    // Find the sample with the greatest distant to the local extreme.
    double max_offset = 0.0;
    for (int i = left_i; i < right_i; ++i)
    {
      double offset = usable_samples[i]->w() - extreme_w;
      if (std::abs(offset) > std::abs(max_offset))
        max_offset = offset;
    }

    // Request an extra sample on the other side of the critical point than that we already have.
    HorizontalDirection direction = max_offset > 0.0 ? HorizontalDirection::left : HorizontalDirection::right;
    w = extreme_w + cubic_used_->scale().step(direction);

    // Remember this direction in hdirection_. This is used when we can't find any local minimum in the fourth degree polynomial.
    if (hdirection_ == HorizontalDirection::undecided)
    {
      hdirection_ = direction;
      Dout(dc::notice, "Set hdirection_ to " << hdirection_ << ".");
      Debug(event_server_.trigger(AlgorithmEventType{hdirection_known_event, *last_extreme_, hdirection_}));
    }

    Dout(dc::notice, "Not enough samples to fit a fourth degree polynomial: asking for another samples at w = " << w);
    state_ = IterationState::extra_sample;
    return true;
  }

  // After finding a maximum we want to find a minimum and visa versa. Change next_extreme_type_.
  next_extreme_type_ = opposite(next_extreme_type_);
  Dout(dc::notice, "next_extreme_type_ is toggled to " << next_extreme_type_ << ".");

  // Use a real Sample, not a fake one.
  Sample const& w2 = *usable_samples[local_extreme_index];
  double const w2_1 = w2.w();

  // Brute force find the two samples that, together with the local extreme sample w2, have the largest spread.
  int i0 = left_i != local_extreme_index ? left_i : left_i + 1;
  int i1 = i0 + 1 != local_extreme_index ? i0 + 1 : i0 + 2;
  double best_spread = 0.0;
  if (number_of_usable_samples > 3)
  {
    for (int t0 = i0; t0 < right_i; ++t0)
    {
      if (t0 == local_extreme_index)
        continue;
      for (int t1 = t0 + 1; t1 < right_i; ++t1)
      {
        if (t1 == local_extreme_index)
          continue;
        double w0_1 = usable_samples[t0]->w();
        double w1_1 = usable_samples[t1]->w();
        double spread = utils::square(w0_1 - w1_1) + utils::square(w0_1 - w2_1) + utils::square(w1_1 - w2_1);
        if (spread > best_spread)
        {
          best_spread = spread;
          i0 = t0;
          i1 = t1;
        }
      }
    }
  }

  Sample const& w1 = *usable_samples[i0];
  Sample const& w0 = *usable_samples[i1];

  double w2_2 = w2_1 * w2_1;
  double w2_3 = w2_2 * w2_1;
  double w2_4 = w2_2 * w2_2;

  double w1_1 = w1.w();
  double w1_2 = w1_1 * w1_1;
  double w1_3 = w1_2 * w1_1;
  double w1_4 = w1_2 * w1_2;

  double w0_1 = w0.w();
  double w0_2 = w0_1 * w0_1;
  double w0_3 = w0_2 * w0_1;
  double w0_4 = w0_2 * w0_2;

  Dout(dc::notice, "Fitting a fourth degree polynomial using the samples at " << w0_1 << ", " << w1_1 << " and " << w2_1);

  // If we have (at least) three points, the approximation is a fourth degree polynomial.
  //
  //   A(w) = a + b w + c w² + d w³ + e w⁴
  //
  // for which we determined the value of the derivative at three points, L'(w₀), L'(w₁) and L'(w₂).
  //
  // The matrix form for the coefficients of the polynomial then becomes:
  //
  //   ⎡  1      2w₂      3w₂²     4w₂³ ⎤ ⎡b⎤   ⎡    L'(w₂)   ⎤
  //   ⎢  1      2w₁      3w₁²     4w₁³ ⎥ ⎢c⎥   ⎢    L'(w₁)   ⎥
  //   ⎢  1      2w₀      3w₀²     4w₀³ ⎥ ⎢d⎥ = ⎢    L'(w₀)   ⎥
  //   ⎣w₂-w₁  w₂²-w₁²  w₂³-w₁³  w₂⁴-w₁⁴⎦ ⎣e⎦   ⎣L(w₂) - L(w₁)⎦
  //
  Eigen::Matrix4d M;
  M <<        1.0,   2.0 *  w2_1,    3.0 * w2_2,    4.0 * w2_3,
              1.0,   2.0 *  w1_1,    3.0 * w1_2,    4.0 * w1_3,
              1.0,   2.0 *  w0_1,    3.0 * w0_2,    4.0 * w0_3,
      w2_1 - w1_1,   w2_2 - w1_2,   w2_3 - w1_3,   w2_4 - w1_4;

  Eigen::Vector4d D;
  D <<         w2.dLdw(),
               w1.dLdw(),
               w0.dLdw(),
       w2.Lw() - w1.Lw();

  Eigen::Vector4d C = M.colPivHouseholderQr().solve(D);

  math::Polynomial fourth_degree_approximation(5 COMMA_CWDEBUG_ONLY("w"));
  fourth_degree_approximation[1] = C[0];
  fourth_degree_approximation[2] = C[1];
  fourth_degree_approximation[3] = C[2];
  fourth_degree_approximation[4] = C[3];
  fourth_degree_approximation[0] = w2.Lw() - fourth_degree_approximation(w2_1);

  Dout(dc::notice, "approximation = " << fourth_degree_approximation);

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{fourth_degree_approximation_event, fourth_degree_approximation});
#endif

  auto derivative = fourth_degree_approximation.derivative();
  Dout(dc::notice, "derivative = " << derivative);

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{derivative_event, derivative});
#endif

  double remainder;
  auto quotient = derivative.long_division(extreme_w, remainder);
  Dout(dc::notice, "quotient = " << quotient << " with remainder " << remainder);

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{quotient_event, quotient});
#endif

  std::array<double, 2> zeroes;
  int number_of_zeroes = quotient.get_roots(zeroes);

  // It is possible that the zeroes are not usable because they are on the wrong side.
  if (hdirection_ != HorizontalDirection::undecided)
  {
    auto wrong_side = [this, w](double zero) { return (hdirection_ == HorizontalDirection::left) != (zero < w); };
    for (int zero = 0; zero < number_of_zeroes;)
      if (wrong_side(zeroes[zero]))
      {
        if (--number_of_zeroes == 1 && zero == 0)
          zeroes[0] = zeroes[1];
      }
      else
        ++zero;
  }

  if (number_of_zeroes > 1)
    Dout(dc::notice, "with other local extremes at " << zeroes[0] << " and " << zeroes[1]);
  else if (number_of_zeroes == 1)
    Dout(dc::notice, "with one other local extreme at " << zeroes[0]);
  else
    Dout(dc::notice, "with no other local extremes!");

  if (number_of_zeroes > 0)
  {
    std::array<double, 2> expected_Lw;
    int best_zero = 0;
    if (number_of_zeroes == 2)
    {
      if (extreme_w > zeroes[0] == extreme_w > zeroes[1])
      {
        // If both zeroes are on the same side of the found local extreme, pick the nearest one.
        best_zero = std::abs(extreme_w - zeroes[0]) < std::abs(extreme_w - zeroes[1]) ? 0 : 1;
      }
      else
      {
        for (int zero = 0; zero < number_of_zeroes; ++zero)
          expected_Lw[zero] = fourth_degree_approximation(zeroes[zero]);
        best_zero = (number_of_zeroes == 2 &&
            (hdirection_ == HorizontalDirection::right ||
             (hdirection_ == HorizontalDirection::undecided &&
              expected_Lw[1] < expected_Lw[0]))) ? 1 : 0;
      }
    }
    else
    {
      ASSERT(number_of_zeroes == 1);    // Huh?
      // There is only one usable zero.
      expected_Lw[0] = fourth_degree_approximation(zeroes[0]);
    }
    w = zeroes[best_zero];
    expected_Lw_ = expected_Lw[best_zero];
    have_expected_Lw_ = true;
//    reset_history();
    left_of_ = chain_.end();
    right_of_ = chain_.end();
#ifdef CWDEBUG
    event_server_.trigger(AlgorithmEventType{left_of_right_of_event,
        nullptr, nullptr, next_extreme_type_, w, cubic_used_->cubic()(w)});
#endif

    Dout(dc::notice(hdirection_ == HorizontalDirection::undecided && number_of_zeroes == 2),
        "Best other local extreme was " << zeroes[best_zero] << " with A(" << zeroes[best_zero] << ") = " <<
        fourth_degree_approximation(zeroes[best_zero]) <<
        " (the other has value " << fourth_degree_approximation(zeroes[1 - best_zero]) << ")");
  }
  else
  {
    // Keep going in the same hdirection.
#ifdef CWDEBUG
    bool hdirection_is_undecided = hdirection_ == HorizontalDirection::undecided;
#endif

    // Note: if hdirection_ is still undecided at this point then this will make a
    // step in the direction from the sample (part of scale) that is the furthest away
    // from the critical point of the approximation, towards that critical point
    // and updates hdirection_ accordingly.
    w = last_extreme_->w() + cubic_used_->scale().step(hdirection_);
    expected_Lw_ = cubic_used_->cubic()(w);
    have_expected_Lw_ = true;

    Dout(dc::notice(hdirection_is_undecided), "Initialized hdirection_ to " << hdirection_ << ".");
  }

  if (hdirection_ == HorizontalDirection::undecided)
  {
    // Now that the decision on which hdirection_ we explore is taken, store that decision.
    hdirection_ = w < w2_1 ? HorizontalDirection::left : HorizontalDirection::right;
    Dout(dc::notice, "Initialized hdirection_ to " << hdirection_ << ".");
  }

  // Move hdirection arrow to new extreme.
  Debug(event_server_.trigger(AlgorithmEventType{hdirection_known_event, *last_extreme_, hdirection_}));

#if 0 //FIXME: add this functionality (copied from histogram.cxx)
  // If this is an extreme that we found by exploring hdirection_ from a previous extreme,
  // then mark that last extreme as being explored in that hdirection_.
  if (std::abs(last_step_) == 1)
  {
    mark_explored(last_w_, hdirection_);

    // If the current extreme is a maximum, then it can be marked as having been explored to the opposite direction (in that case the previous extreme is a minimum).
    // Or if the current extreme is a minimum and we came from another minimum (with steps of size 1), then mark it as having been explored into the opposite direction.
    if (saw_minimum_)
    {
      // Mark that the current extreme is being explored in the opposite hdirection_.
      mark_explored(w, opposite(hdirection_));
    }
  }
#endif

  // Keep track of the best minimum so far; or abort if this minimum isn't better than one found before.
  if (last_extreme_->get_extreme_type() == ExtremeType::minimum)
  {
    //saw_minimum_ = true;
    if (best_minimum_ == chain_.end() || best_minimum_->Lw() > last_extreme_->Lw())
      best_minimum_ = last_extreme_;
    if (last_extreme_ != best_minimum_)
      return false;
  }

  state_ = IterationState::find_extreme;
  return true;
}

} // namespace gradient_descent
