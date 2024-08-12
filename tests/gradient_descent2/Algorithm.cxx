#include "sys.h"
#include "AnalyzedCubic.h"
#include "SampleNode.h"
#include "Algorithm.h"
#include "IterationState.h"
#include "debug.h"
#ifdef CWDEBUG
#include "utils/at_scope_end.h"
#endif

namespace gradient_descent {

bool Algorithm::operator()(double& w, double Lw, double dLdw)
{
  DoutEntering(dc::notice, "Algorithm::operator()(" << w << ", " << Lw << ", " << dLdw << ")");

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

  // Print the chain_ to debug output just before leaving this function.
  auto&& dump_chain = at_scope_end([this]{ chain_.dump(); });
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
      Dout(dc::notice, "Running " << state_);
      switch (state_)
      {
        case IterationState::first_cubic:
        {
          auto left_node = new_node;
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
          left_node->initialize_cubic(right_node COMMA_CWDEBUG_ONLY(event_server_, current_is_left));

          // Jump to extreme of the first cubic.
          w = left_node->find_extreme(*right_node, next_extreme_type_);
          cubic_used_ = left_node;

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

            cubic_used_->set_scale(scale_cp_type, critical_point_w, left_node->w(), right_node->w());
#ifdef CWDEBUG
            event_server_.trigger(AlgorithmEventType{scale_draw_event, ScaleUpdate::initialized, *left_node, {}});
#endif
          }

          // If cubic_used_ didn't have a usable extreme, then next_extreme_type_ was set to unknown.
          if (next_extreme_type_ == ExtremeType::unknown)
          {
            // We're going to look for a minimum.
            if (cubic_used_->is_rising())
            {
              // Keep going to the left.
              w = left_node->w() - cubic_used_->scale().value();
              left_of_ = left_node;
            }
            else // cubic_used_->is_falling() or the cubic is flat.
            {
              // Keep going to the right.
              w = right_node->w() + cubic_used_->scale().value();
              if (cubic_used_->is_falling())
                right_of_ = right_node;
            }
          }
          else
          {
            ASSERT(w != SampleNode::uninitialized_magic);
            double const negligible_offset = negligible_scale_fraction * cubic_used_->scale().value();
            if (left_node->w() - w > negligible_offset)
              left_of_ = left_node;
            else if (right_node->w() - w > negligible_offset)
            {
              left_of_ = right_node;
              if (w - left_node->w() > negligible_offset)
                right_of_ = left_node;
            }
            else if (w - right_node->w() > negligible_offset)
              right_of_ = right_node;
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
          if (std::abs(w - new_node->w()) <= 1e-3 * left_node->scale().value())     // Too close to the last sample?
            if (!handle_local_extreme(w))
              return false;

          break;
        }

        case IterationState::find_extreme:
        {
#ifdef CWDEBUG
          math::CubicPolynomial old_cubic = cubic_used_ != chain_.end() ? cubic_used_->cubic() : math::CubicPolynomial{};
#endif
          auto right_node = std::next(new_node);

          AnalyzedCubic right_acubic;       // If right_node != chain_.end() is true then this corresponds to the cubic of std::next(new_node).
          bool right_has_extreme = false;
          double right_extreme_w;
          double right_dist;
          if (right_node != chain_.end())
          {
            new_node->initialize_cubic(right_node COMMA_CWDEBUG_ONLY(event_server_, true));        // true: called on new_node.
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
            auto left_node = std::prev(new_node);
            left_node->initialize_cubic(new_node COMMA_CWDEBUG_ONLY(event_server_, false));        // false: new_node is passed.
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
          if (step < (next_extreme_type_ == ExtremeType::minimum ? 0.01 : 0.05) * cubic_used_->scale().value())
            if (!handle_local_extreme(w))
              return false;

          break;
        }

        case IterationState::finish:
        {
          // Store the last sample and set the state to success.
          set_global_minimum(new_node);
          return false;
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

  // FIXME: implement.
  ASSERT(false);
  return false;
}

bool Algorithm::handle_local_extreme(double& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_local_extreme(" << w << ")");

  // Set chain_.larger_ to point to the first SampleNode that is larger than w, if any.
  chain_.find_larger(w);

  // Find out if a new (fake) SampleNode should inserted.
  SampleNode::const_iterator right = chain_.larger();
  SampleNode::const_iterator left = right == chain_.begin() ? chain_.end() : std::prev(right);

  double const negligible_offset = negligible_scale_fraction * cubic_used_->scale().value();

  SampleNode::const_iterator local_extreme;
  if (right != chain_.end() && std::abs(right->w() - w) <= negligible_offset)
  {
    Dout(dc::notice, "Using existing sample on the right at " << right->w());
    local_extreme = right;
  }
  else if (left != chain_.end() && std::abs(left->w() - w) <= negligible_offset)
  {
    Dout(dc::notice, "Using existing sample on the left at " << left->w());
    local_extreme = left;
  }
  else
  {
    // Create a new Sample for this sample.
    Sample fake_sample{w, cubic_used_->cubic()(w), 0.0 COMMA_CWDEBUG_ONLY(label_context_)};
    Dout(dc::notice, "Using fake sample = " << fake_sample);
    local_extreme = chain_.insert(std::move(fake_sample));
    local_extreme->set_fake(true);
    auto next = std::next(local_extreme);
    if (next != chain_.end())
      local_extreme->initialize_cubic(next, event_server_, true);
    if (local_extreme != chain_.begin())
      std::prev(local_extreme)->initialize_cubic(local_extreme COMMA_CWDEBUG_ONLY(event_server_, false));
  }

#ifdef CWDEBUG
  // Draw the new sample.
  event_server_.trigger(AlgorithmEventType{new_local_extreme_event, *local_extreme, std::to_string(local_extreme->label())});
#endif

  finish();
  return true;
}

} // namespace gradient_descent
