#include "sys.h"
#include "AnalyzedCubic.h"
#include "SampleNode.h"
#include "Algorithm.h"
#include "IterationState.h"
#include <Eigen/Dense>
#include "debug.h"
#ifdef CWDEBUG
#include "utils/at_scope_end.h"
#include "utils/print_using.h"
#include "cwds/Restart.h"
#endif

namespace gradient_descent {

void Algorithm::initialize_node(SampleNode::iterator node, SampleNode::const_iterator next
    COMMA_CWDEBUG_ONLY(bool node_is_last))
{
  // Remember if node (the left node of the cubic) is already marked as an extreme.
  // If that is the case then it will have either right_max_bit or right_min_bit set.
  bool const is_extreme_before = (static_cast<int>(node->type()) & (right_max_bit|right_min_bit)) != 0;

  // Initialize the cubic between node and next.
  node->initialize_cubic(next
      // If state_ is extra_sample, then pass ExtremeType::unknown just to avoid an assert.
      COMMA_CWDEBUG_ONLY(state_ == IterationState::extra_sample ? ExtremeType::unknown : next_extreme_type_, event_server_, node_is_last));

  bool const is_extreme_after = (static_cast<int>(node->type()) & (right_max_bit|right_min_bit)) != 0;

  // If node was already marked as extreme, then it should still be!
  ASSERT(!is_extreme_before || is_extreme_after);

  if (!is_extreme_before && is_extreme_after)
  {
    bool const is_maximum_after = (static_cast<int>(node->type()) & right_max_bit) != 0;
    // The node has been marked as right_* extreme type.
    // Mark the opposite side of the extreme as left_* extreme type.
    if (node != chain_.begin())
      std::prev(node)->change_type_to_left_extreme(is_maximum_after ? ExtremeType::maximum : ExtremeType::minimum);
  }
}

void Algorithm::fix_flat(SampleNode::iterator left_node, SampleNode::iterator fix_node)
{
  // Assume that fix_node is a minimum:
  //
  //                left_node  right_node
  //                     |         |
  //                    \v         v/
  //                     O         O
  //  left_node-fix_node →\       /← fix_node-right_node
  //         cubic         \     /         cubic
  //                        `-O-´
  //                      __↗ ^ ↖__
  //                     /    |    \
  //                    /  fix_node \
  //                   /             (*)
  // Is the right part of the "left_node-fix_node cubic" flat?
  CubicEndShape left_side_shape = get_end(left_node->type(), false/*right side*/);
  //                                   (*)
  //                                     \
  // Is the left part of the "fix_node-right_node cubic" flat?
  CubicEndShape right_side_shape = get_end(fix_node->type(), true/*left side*/);

  bool left_side_is_flat = left_side_shape == flat_high || left_side_shape == flat_low;
  bool right_side_is_flat = right_side_shape == flat_high || right_side_shape == flat_low;

  if (left_side_is_flat != right_side_is_flat)
  {
    if (left_side_is_flat)
      // Make sure the left part of the "fix_node-right_node cubic" is flat too.
      // Note that change_type_to_right_extreme changes the type into "the right part of an extreme", and
      // thus changes the left end of the type of the cubic.
      fix_node->change_type_to_right_extreme(left_side_shape == flat_high ? ExtremeType::maximum : ExtremeType::minimum);
    else
      // Make sure the right part of the "left_node-fix_node cubic" is flat too.
      left_node->change_type_to_left_extreme(right_side_shape == flat_high ? ExtremeType::maximum : ExtremeType::minimum);
  }
}

bool Algorithm::operator()(WeightRef w, double Lw, double dLdw)
{
  DoutEntering(dc::notice, "Algorithm::operator()(" << w << ", " << Lw << ", " << dLdw << ")");
  RESTART

  // Do not re-enter this function after it returned false!
  ASSERT(state_ != IterationState::success);

  // If the function is uncalculatable.
  bool must_abort = !std::isfinite(Lw) || !std::isfinite(dLdw);

#ifdef CWDEBUG
  // Erase all previous curves (if they exist).
  event_server_.trigger(AlgorithmEventType{reset_event});
  if (!must_abort && have_expected_Lw_)
  {
    // Draw an arrow from where expected to end up, based on the approximation, and the actual value of L(w).
    event_server_.trigger(AlgorithmEventType{difference_event, w, expected_Lw_, Lw});
    have_expected_Lw_ = false;
  }

  // Print the chain_ to debug output and do a sanity check just before leaving this function.
  auto&& dump_chain = at_scope_end([this]{
      // After the last step the sanity check might fail because we skip initializing the new cubics left and right of the final minimum.
      if (state_ != IterationState::success)
      {
        chain_.dump(this);
        chain_.sanity_check(this);
      }
  });
#endif

  SampleNode::const_iterator new_node;
  SampleNode::const_iterator left_node;         // The node left of the new node (only valid if new_node != chain_.begin()).
  SampleNode::const_iterator right_node;        // The node right of the new node.

  if (!must_abort)
  {
    SampleNode::iterator non_const_left_node;
    SampleNode::iterator non_const_new_node;
    SampleNode::iterator non_const_right_node;

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
    non_const_new_node = chain_.insert(std::move(current));
    new_node = non_const_new_node;

    // There is no NEED to initialize the cubics left and right of the final sample.
    // In fact, that would cause an assertion, so avoid it.
    if (AI_UNLIKELY(state_ != IterationState::finish))
    {
      non_const_right_node = std::next(non_const_new_node);
      right_node = non_const_right_node;

      // Initialize the cubic(s).
      if (new_node != chain_.begin())
      {
        non_const_left_node = std::prev(non_const_new_node);
        initialize_node(non_const_left_node, new_node COMMA_CWDEBUG_ONLY(false));         // false: new_node is passed as second argument.
        left_node = non_const_left_node;

        if (non_const_left_node != chain_.begin())
          fix_flat(std::prev(non_const_left_node), non_const_left_node);
      }
      if (right_node != chain_.end())
      {
        initialize_node(non_const_new_node, right_node COMMA_CWDEBUG_ONLY(true));         // true: new_node is passed as first argument.
        if (new_node != chain_.begin()) // Otherwise non_const_left_node wasn't even initialized.
          fix_flat(non_const_left_node, non_const_new_node);
        if (non_const_right_node->is_cubic())
          fix_flat(non_const_new_node, non_const_right_node);
      }
    }
  }

  for (;;)
  {
    // Update kinetic energy. Returns false if too much energy was used.
    if (must_abort || !update_energy(Lw))
    {
      if (!handle_abort_hdirection(w))
        return false;
    }
    else
    {
      Dout(dc::notice, "⁂  Running \"" << state_ << '"');
      RESTART
      switch (state_)
      {
        case IterationState::first_cubic:
        {
          cubic_used_ = right_node == chain_.end() ? left_node : new_node;
          auto next = std::next(cubic_used_);

          // Jump to extreme of the first cubic. This sets next_extreme_type_ to the found/best extreme type, if any.
          w.find_extreme_jump(cubic_used_, *next, next_extreme_type_);

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
              left_of_ = cubic_used_;
              // Keep going to the left; jump to left_of_ and step further to the left by one scale value.
              w.node_jump_step(left_of_, -cubic_used_->scale().value());
            }
            else // cubic_used_->is_falling() or the cubic is flat.
            {
              // Keep going to the right; jump to next and step further to the right by one scale value.
              w.node_jump_step(next, cubic_used_->scale().value());
              if (cubic_used_->is_falling())
                right_of_ = next;
            }
            // Abort this hdirection if we used too much energy due to this step.
            check_energy_ = hdirection_ != HorizontalDirection::undecided;
            Dout(dc::notice(check_energy_), "Set check_energy_ to true because " << *cubic_used_ <<
                " doesn't have a usable extreme and we just added the scale.");
          }
          else
          {
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
              right_of_ == chain_.end() ? nullptr : &*right_of_});
          event_server_.trigger(AlgorithmEventType{jump_point_event,
              next_extreme_type_, w, cubic_used_->cubic()(w)});
#endif

          // If next_extreme_type_ was not set then this cubic has no extrema. In that case go for a minimum.
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
            if (right_acubic.has_extrema())
            {
              //FIXME: initialize scale immediately after an initialize_cubic? And then use the real scale here.
              double const significant_offset = significant_scale_fraction * new_node->w_scale_estimate();
              right_extreme_w = right_acubic.get_extreme();
              right_dist = std::abs(right_extreme_w - new_node->w()) + std::abs(right_extreme_w - right_node->w());
              double const humongous_offset = humongous_step_scale_factor * new_node->w_scale_estimate();
              right_has_extreme =
                right_dist < humongous_offset &&
                (left_of_ == chain_.end() || right_extreme_w - left_of_->w() < significant_offset) &&
                (right_of_ == chain_.end() || right_of_->w() - right_extreme_w < significant_offset);
            }
          }
          AnalyzedCubic left_acubic;        // If new_node != chain_.begin() is true then this corresponds to the cubic of std::prev(new_node).
          bool left_has_extreme = false;
          double left_extreme_w;
          double left_dist;
          if (new_node != chain_.begin())
          {
            left_acubic.initialize(left_node->cubic(), next_extreme_type_);
            if (left_acubic.has_extrema())
            {
              //FIXME: initialize scale immediately after an initialize_cubic? And then use the real scale here.
              double const significant_offset = significant_scale_fraction * left_node->w_scale_estimate();
              left_extreme_w = left_acubic.get_extreme();
              left_dist = std::abs(left_extreme_w - left_node->w()) + std::abs(left_extreme_w - new_node->w());
              double const humongous_offset = humongous_step_scale_factor * left_node->w_scale_estimate();
              left_has_extreme =
                left_dist < humongous_offset &&
                (left_of_ == chain_.end() || left_extreme_w - left_of_->w() < significant_offset) &&
                (right_of_ == chain_.end() || right_of_->w() - left_extreme_w < significant_offset);
            }
          }

          cubic_used_ = chain_.end();
          AnalyzedCubic* used_acubic = nullptr;
          SampleNode::const_iterator used_cubic;        // Only valid if used_acubic was set.

          // last_extreme_cubic_ and hdirection_ are set at the same time (in handle_local_extreme).
          // In the code below, where we test 'hdirection_ == HorizontalDirection::undecided' we assume
          // that if that is false then also last_extreme_cubic_ will be valid.
          ASSERT((last_extreme_cubic_ == chain_.end()) == (hdirection_ == HorizontalDirection::undecided));

          if (right_has_extreme)
          {
            // If both cubics have the required extreme in the required region, then
            // use the cubic that is closest to the last local extreme.
            // If hdirection is still undecided there is no last local extreme yet; in that case use the cubic
            // whose samples lay closest to their extreme.
            if (!left_has_extreme || hdirection_ == HorizontalDirection::left ||
                (hdirection_ == HorizontalDirection::undecided && right_dist < left_dist))
            {
#ifdef CWDEBUG
              Dout(dc::notice|continued_cf, "Choosing " << next_extreme_type_ << " of right cubic because ");
              if (left_has_extreme)
              {
                if (hdirection_ == HorizontalDirection::undecided)
                  Dout(dc::finish, "that extreme (" << right_extreme_w << ") is closer to the center of [" <<
                      new_node->label() << "]<--->[" << right_node->label() << "], than the left cubic extreme (" <<
                      left_extreme_w << ") is from the center of [" << std::prev(new_node)->label() << "]<--->[" << new_node->label() <<
                      "]; (" << right_dist << " < " << left_dist << ").");
                else
                  Dout(dc::finish, "that cubic [" << new_node->label() << "] is closer to the last extreme [" <<
                      last_extreme_cubic_->label() << "] than the left cubic [" << left_node->label() << "].");
              }
              else if (new_node != chain_.begin())
                Dout(dc::finish, "the left one does not have that extreme between " <<
                    (right_of_ != chain_.end() ? right_of_->w() : -std::numeric_limits<double>::infinity()) << " and " <<
                    (left_of_ != chain_.end() ? left_of_->w() : std::numeric_limits<double>::infinity()));
              else
                Dout(dc::finish, "the left one does not exists.");
#endif
              used_cubic = new_node;
              used_acubic = &right_acubic;
            }
            else
            {
#ifdef CWDEBUG
              Dout(dc::notice|continued_cf, "Choosing " << next_extreme_type_ << " of left cubic because ");
              if (left_has_extreme)
              {
                if (hdirection_ == HorizontalDirection::undecided)
                  Dout(dc::finish, "that extreme (" << left_extreme_w << ") is closer to the center of [" <<
                      std::prev(new_node)->label() << "]<--->[" << new_node->label() << "], than the right cubic extreme (" <<
                      right_extreme_w << ") is from the center of [" << new_node->label() << "]<--->[" << right_node->label() <<
                      "]; (" << left_dist << " < " << right_dist << ").");
                else
                  Dout(dc::finish, "that cubic [" << left_node->label() << "] is closer to the last extreme [" <<
                      last_extreme_cubic_->label() << "] than the right cubic [" << new_node->label() << "].");
              }
              else
                Dout(dc::finish, "the right one does not have that extreme between " <<
                    (right_of_ != chain_.end() ? right_of_->w() : -std::numeric_limits<double>::infinity()) << " and " <<
                    (left_of_ != chain_.end() ? left_of_->w() : std::numeric_limits<double>::infinity()));
#endif
              used_cubic = std::prev(new_node);
              used_acubic = &left_acubic;
            }
          }
          else if (left_has_extreme)
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
            used_cubic = std::prev(new_node);
            used_acubic = &left_acubic;
          }

          // If used_cubic isn't set then keep_going_node and keep_going_direction will be set.
          SampleNode::const_iterator keep_going_node;                                   // The node that we add scale to get the new jump point.
          HorizontalDirection keep_going_direction = HorizontalDirection::undecided;    // Default: undecided means nothing needs to be done.

          if (left_has_extreme || right_has_extreme)
            handle_extreme_jump(w, used_cubic, *used_acubic, {});
          else
          {
#ifdef CWDEBUG
            if (new_node == chain_.begin())
            {
              ASSERT(right_node != chain_.end());
              Dout(dc::notice|continued_cf, "There is no left cubic, and the right cubic has no ");
            }
            else if (right_node == chain_.end())
              Dout(dc::notice|continued_cf, "There is no right cubic, and the left cubic has no ");
            else
              Dout(dc::notice|continued_cf, "Neither the left nor the right cubic has a ");
            Dout(dc::finish, next_extreme_type_ << " in the range " << utils::print_using(*this, &Algorithm::print_range_on) << '.');
#endif
            SampleNode::const_iterator current = new_node;
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
            // less than zero due to floating-point round-off errors).

            HorizontalDirection direction;
            if (left_of_ == chain_.end() && right_of_ == chain_.end())
            {
              // If no region boundaries are set, then the cubics don't have an extreme at all; aka,
              // we are in the middle of a monotonic rising or falling area: in that case just
              // move the current point to the highest/lowest point left/right of the current sample,
              // depending on whether we're looking for a minimum or maximum.
              bool at_right_edge = !current->is_cubic();
              if (at_right_edge)
                current = std::prev(current);
              ASSERT(current->is_cubic());
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
              while (current->is_cubic() && !current->has_unfound_extreme(next_extreme_type_))
                current = std::next(current);
              if (!current->is_cubic())
              {
                right_of_ = current;
                // Set cubic_used_ in order to set the scale on it: if we keep going again, we want to use the scale of this cubic.
                cubic_used_ = std::prev(current);
                right_acubic.initialize(cubic_used_->cubic(), next_extreme_type_);
                used_acubic = &right_acubic;
                // Keep going to the right. cubic_used_->scale().value() needs to be added to right_of_->w().
                keep_going_node = right_of_;
                keep_going_direction = HorizontalDirection::right;
              }
              else
              {
                // implement
                ASSERT(false);
              }
            }
            else
            {
              while (current != chain_.begin() && !std::prev(current)->has_unfound_extreme(next_extreme_type_))
                current = std::prev(current);
              if (current == chain_.begin())
              {
                left_of_ = current;
                // Set cubic_used_ in order to set the scale on it: if we keep going again,
                // we want to use the scale of this cubic.
                cubic_used_ = current;
                left_acubic.initialize(cubic_used_->cubic(), next_extreme_type_);
                used_acubic = &left_acubic;
                // Keep going to the left. cubic_used_->scale().value() needs to be subtracted from left_of_->w().
                keep_going_node = left_of_;
                keep_going_direction = HorizontalDirection::left;
              }
              else
              {
                // implement
                ASSERT(false);
              }
            }
          }

          // cubic_used_ should have been set to the cubic that we need to update the scale of.
          ASSERT(cubic_used_ != chain_.end());

          // Update scale.
          CriticalPointType const scale_cp_type =
            used_acubic->has_extrema() ? (next_extreme_type_ == ExtremeType::minimum ? CriticalPointType::minimum
                                                                                      : CriticalPointType::maximum)
                                       : CriticalPointType::inflection_point;

          double const critical_point_w =
            used_acubic->has_extrema() ? used_acubic->get_extreme()
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
          event_server_.trigger(AlgorithmEventType{scale_draw_event, ScaleUpdate::away_from_cp, *cubic_used_, old_cubic});
#endif

          if (keep_going_direction != HorizontalDirection::undecided)
          {
            // Now that scale is known, add scale if we must "keep going" (there was no extreme to jump to).
            w.node_jump_step(keep_going_node, static_cast<int>(keep_going_direction) * cubic_used_->scale().value());
          }
          else
          {
            // Adjust left_of_/right_of_.
            double const negligible_offset = negligible_scale_fraction * cubic_used_->scale().value();
            if (new_node->w() - w > negligible_offset)
              left_of_ = new_node;
            else if (w - new_node->w() > negligible_offset)
              right_of_ = new_node;
          }

          //FIXME: what to do in the case where this fails?
          ASSERT(debug_within_range(w));

#ifdef CWDEBUG
          event_server_.trigger(AlgorithmEventType{left_of_right_of_event,
              left_of_ == chain_.end() ? nullptr : &*left_of_,
              right_of_ == chain_.end() ? nullptr : &*right_of_});
          event_server_.trigger(AlgorithmEventType{jump_point_event,
              next_extreme_type_, w, cubic_used_->cubic()(w)});
#endif

          // If the last step is significantly smaller than the scale, then we found a local extreme.
          double step = std::abs(chain_.last()->w() - w);
          Dout(dc::notice, "step = " << step);
          if (step < (next_extreme_type_ == ExtremeType::minimum ? 0.01 : 0.05) * cubic_used_->scale().value())
          {
            if (!handle_local_extreme(w) && !handle_abort_hdirection(w))
              return false;
          }
          else
          {
            // Abort this hdirection if we used too much energy due to this step.
            check_energy_ = keep_going_direction != HorizontalDirection::undecided && hdirection_ != HorizontalDirection::undecided;
            Dout(dc::notice(check_energy_), "Set check_energy_ to true because " <<
                *cubic_used_ << " had no extreme to jump to and we just added the scale.");
          }

          break;
        }

        case IterationState::extra_sample:
        {
          if (!handle_local_extreme(w) && !handle_abort_hdirection(w))
            return false;

          break;
        }

        case IterationState::backtracking:
        {
          // If the new sample results in a cubic, adjacent to the last local extreme, that shouldn't be used,
          // then keep_going_node and keep_going_direction will be set.
          SampleNode::const_iterator keep_going_node;                                   // The node that we add scale to get the new jump point.
          HorizontalDirection keep_going_direction = HorizontalDirection::undecided;    // Default: undecided means nothing needs to be done.

          // The left- or right-most sample used to fit the fourth degree polynomial for the last local extremum found.
          SampleNode::const_iterator local_extreme_edge;        // Only valid if we have a last_extreme_cubic_.

          if (hdirection_ != HorizontalDirection::undecided)
          {
            // Initialize local_extreme_edge.
            local_extreme_edge = last_extreme_cubic_->local_extreme().get_edge_node(hdirection_);
            bool const adjacent_to_last_extreme =
              ((hdirection_ == HorizontalDirection::left && local_extreme_edge == right_node) ||
              (hdirection_ == HorizontalDirection::right && local_extreme_edge == left_node));

            if (adjacent_to_last_extreme)
            {
              SampleNode::const_iterator adjacent_cubic = hdirection_ == HorizontalDirection::left ? new_node : left_node;
              AnalyzedCubic adjacent_acubic;
              adjacent_acubic.initialize(adjacent_cubic->cubic(), opposite(next_extreme_type_));

              // The expected end of the adjacent_cubic if it is heading towards the sought local extremum.
              CubicEndShape const correct_end =
                (hdirection_ == HorizontalDirection::right) == (next_extreme_type_ == ExtremeType::minimum) ? downhill : uphill;

              // The cubic ends with a slope as expected if heading towards the next sought-for local extremum.
              bool const has_expected_slope_for_approaching_sought_for_local_extreme =
                get_end(adjacent_cubic->type(), hdirection_ == HorizontalDirection::left) == correct_end;

              bool const is_approaching_sought_for_local_extreme =
                // Does not contain the next sought local extremum.
                !adjacent_cubic->has_unfound_extreme(next_extreme_type_) &&
                has_expected_slope_for_approaching_sought_for_local_extreme;

              bool const is_approaching_the_wrong_local_extreme =
                !adjacent_cubic->has_unfound_extreme(next_extreme_type_) &&
                !has_expected_slope_for_approaching_sought_for_local_extreme;

              // Get a coordinate that is further away than where the real local extreme of the last local extreme could be.
              double const significantly_far_away = local_extreme_edge->w() + static_cast<int>(hdirection_) *
                (next_extreme_type_ == ExtremeType::minimum ? 0.055 : 0.011) * last_extreme_cubic_->scale().value();

              // The cubic first seems to pass an extreme of the opposite type.
              bool const definitely_skipped_next_local_extreme =
                is_approaching_sought_for_local_extreme &&
                adjacent_cubic->has_extreme(opposite(next_extreme_type_)) &&
                static_cast<int>(hdirection_) * (adjacent_acubic.get_extreme() - significantly_far_away) > 0.0;

              if (definitely_skipped_next_local_extreme || is_approaching_the_wrong_local_extreme)
              {
#ifdef CWDEBUG
                if (definitely_skipped_next_local_extreme)
                  Dout(dc::notice, "Not using this cubic however, because it contains again and only a " << opposite(next_extreme_type_) <<
                      ", meaning it completely missed the next local extreme.");
                else
                  Dout(dc::notice, "Ignoring cubic adjacent to last local extreme (" << last_extreme_cubic_->local_extreme().label() <<
                      " of type " << last_extreme_cubic_->local_extreme().get_extreme_type() << ") because it is heading for a " <<
                      opposite(next_extreme_type_) << " again!");
#endif
                keep_going_direction = hdirection_;
              }
              else
              {
                // So that we can use break;
                for (int once = 0; once != 1; once = 1)
                {
                  // Let (x1, y1) be the point where the approximation switches from the
                  // fourth degree polynomial to the new cubic (local_extreme_edge).
                  // Let P(x) be the fourth degree polynomial and let Q(x) be the cubic (adjacent_cubic).
                  // Thus P(x1) = Q(x1) and P'(x1) = Q'(x1) (the piece-wise approximation is C1 continuous).
                  //
                  double x1 = local_extreme_edge->w();
                  double y1 = local_extreme_edge->Lw();
                  // Let x2 be a point a scale distance into Q.
                  double const scale_step = last_extreme_cubic_->scale().step(hdirection_);
                  double x2 = x1 + scale_step;
                  // Break if we stepped over Q instead of into it.
                  if ((hdirection_ == HorizontalDirection::left && x2 < adjacent_cubic->w()) ||
                      (hdirection_ == HorizontalDirection::right && x2 > adjacent_cubic->next_node()->w()))
                    break;
                  math::Polynomial const& fourth_degree_approximation = last_extreme_cubic_->local_extreme().get_fourth_degree_approximation();
                  // Let Py = P(x2)
                  double Py = fourth_degree_approximation(x2);
                  // Calculate Qy = Q(x2).
                  double Qy = adjacent_cubic->cubic()(x2);

                  // Apply a coordinate transformation such that (x1, y1) --> (0, 0) and (x2, Py) --> (1, 1):
                  // Note that x2' = 1, and Py' = 1.

                  // Then Qy' = (Qy - y1) / (Py - y1);
                  double Qy_prime = (Qy - y1) / (Py - y1);

                  // Take dot product between the vectors from (x1, y1) to (x2, Py) and (x2, Qy) respectively, using the transformed coordinates.
                  // Thus x2' * x2' + Py' * Qy' = 1 + Qy'.
                  // The angle between the two vectors then is given by:
                  // 1 + Qy' = |(x2', Py')| |x2', Qy'| cos(angle) = sqrt(2) * sqrt(1 + Qy'^2) cos(angle)
                  // Calculate the square of the cosine of the angle:
                  double two_cos_squared = utils::square(1.0 + Qy_prime) / (1.0 + utils::square(Qy_prime));

                  if (two_cos_squared < 1.8)
                  {
                    Dout(dc::notice, "Not using this cubic however, because the cubic deviates too much from the fourth degree polynomial "
                        "(P(" << x2 << ") = " << Py << " and Q(" << x2 << ") = " << Qy << " --> two_cos_squared = " << two_cos_squared << ".");

                    // This cubic "fit" is likely nonsense. Lets just make a scale step.
                    keep_going_direction = hdirection_;
#ifdef CWDEBUG
                    event_server_.trigger(AlgorithmEventType{fourth_degree_approximation_event, fourth_degree_approximation});
#endif
                  }
#ifdef CWDEBUG
                  else
                    Dout(dc::notice, "two_cos_squared = " << two_cos_squared);
#endif
                }
              }
            }
          }

          // Do not run backtracking again.
          state_ = IterationState::find_extreme;

          if (keep_going_direction == HorizontalDirection::undecided)   // Ok to use the find_extreme algorithm?
            continue;

          keep_going_node = hdirection_ == HorizontalDirection::left ? left_of_ : right_of_;
          ASSERT(keep_going_node != chain_.end());
          cubic_used_ = last_extreme_cubic_;
          ASSERT(cubic_used_ != chain_.end());
          chain_.reuse(local_extreme_edge);
          w.node_jump_step(keep_going_node, static_cast<int>(keep_going_direction) * cubic_used_->scale().value());

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
    auto ibp = chain_.duplicate(cubic_used_->scale().value(), state_ == IterationState::finish);
    if (AI_LIKELY(!ibp.second))
    {
      // Ask for a new sample at `w` if not already finished.
      Dout(dc::notice, "Next probe: " << std::setprecision(std::numeric_limits<long double>::max_digits10) << w);
      break;
    }

    Dout(dc::notice, "New probe (" << w << ") is too close to an existing sample (" << *ibp.first << "), reusing that.");
    new_node = ibp.first;
    chain_.reuse(new_node);

    // Re-initialize right_node and left_node.
    right_node = std::next(new_node);          // The node right of the new node.
    if (new_node != chain_.begin())            // If this is true then left_node isn't used (and therefore doesn't need to be assigned).
      left_node = std::prev(new_node);

    must_abort = false;
  }

  return true;
}

void Algorithm::handle_single_sample(WeightRef w)
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

  Debug(set_algorithm_str(w, algorithm_str));
  handle_single_sample_step(w, step, {});
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

void Algorithm::move_into_range(WeightRef w)
{
  // Before this function is called, we expect that initialize_range was called,
  // which should have searched from the last extreme in the direction of hdirection_
  // until it found a single cubic containing next_extreme_type_.
  //
  // Hence, if such a cubic was found, left_of_ should be equal to right_of_->next_node().
  // Otherwise, if hdirection_ is left, left_of_ should be equal to chain_.begin() and right_of_ equal to chain_.end(),
  // or if hdirection_ is right, left_of_ should be equal to chain_.end() and right_of_ equal to std::prev(chain_.end()).
#if CW_DEBUG
  if (left_of_ != chain_.end() && right_of_ != chain_.end())
    ASSERT(right_of_->next_node() == left_of_);
  else if (hdirection_ == HorizontalDirection::left)
    ASSERT(left_of_ == chain_.begin() && right_of_ == chain_.end());
  else
    ASSERT(hdirection_ == HorizontalDirection::right &&
           left_of_ == chain_.end() && right_of_ == std::prev(chain_.end()));
  // Since hdirection_ is not allowed to be HorizontalDirection::undecided here (already enforced by calling initialize_range)
  // the above assures that at least one of left_of_ and right_of_ is unequal chain_.end().
#endif

  if (left_of_ != chain_.end() && right_of_ != chain_.end())
  {
    // We have a cubic: right_of_ is that cubic:
    //
    // right_of_ (in this case, the left-side of the range is a maximum)
    //    |
    //    v
    //    -.  left_of_ (in this case, the right-side of the range is a minimum)
    //      \   |
    //       \  |
    //        \ v
    //         `-
    //    |<--->|
    //     range

    // Due to floating-point round off errors, it is safer to jump to left_of_/right_of_ respectively if those are flat with the correct extreme.
    if (right_of_->has_flat_extreme(next_extreme_type_, on_the_right_side))
      extreme_jump_range_edge(w, left_of_, {});
    else if (right_of_->has_flat_extreme(next_extreme_type_, on_the_left_side))
      extreme_jump_range_edge(w, right_of_, {});
    // Jump to the extreme of the cubic if w is out of range. Otherwise keep w (do nothing).
    else if (left_of_->w() < w || w < right_of_->w())
    {
      // w is out of range.
      Dout(dc::notice, "Attempt to jump to " << std::setprecision(17) << w <<
          " thwarted because that is not in the required range [" << right_of_->w() << ", " << left_of_->w() << "].");
      AnalyzedCubic acubic;
      acubic.initialize(right_of_->cubic(), next_extreme_type_);
      handle_extreme_jump(w, right_of_, acubic, {});
    }
  }
  else
  {
    bool const too_far_to_the_right = left_of_ != chain_.end() && left_of_->w() < w;
    bool const too_far_to_the_left = right_of_ != chain_.end() && w < right_of_->w();

    // If we're already in range - do nothing.
    if (!too_far_to_the_right && !too_far_to_the_left)
      return;

#ifdef CWDEBUG
    if (too_far_to_the_right)
    {
      Dout(dc::notice, "Attempt to jump to " << w << " thwarted because that is larger than left_of_ (" << left_of_->w() << ").");
      // If only left_of_ is set then we should be going to the left and left_of_ must be chain_.begin().
      ASSERT(hdirection_ == HorizontalDirection::left && left_of_ == chain_.begin());
    }
    else
    {
      Dout(dc::notice, "Attempt to jump to " << w << " thwarted because this is less than right_of_ (" << right_of_->w() << ".");
      // If only right_of_ is set then we should be going to the right and right_of_ shouldn't have a cubic defined.
      ASSERT(hdirection_ == HorizontalDirection::right && !right_of_->is_cubic());
    }
#endif

    // We're out of range and need to make a step in hdirection_.
    Scale const& scale = too_far_to_the_right ? left_of_->scale() : std::prev(right_of_)->scale();
    double step = scale.is_valid() ? scale.step(hdirection_) : static_cast<int>(hdirection_) * small_step_;
    w.node_jump_step(too_far_to_the_right ? left_of_ : right_of_, step);
  }

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{jump_point_event,
      next_extreme_type_, w, expected_Lw_});
#endif

  ASSERT(debug_within_range(w));
}

bool Algorithm::update_energy(double Lw)
{
  if (check_energy_)
  {
    check_energy_ = false;

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
      if (last_extreme_cubic_ != chain_.end())
      {
        // We tried to go to hdirection_ but that wasn't possible.
        // Mark the local extreme that we came from as explored, in this direction.
        last_extreme_cubic_->local_extreme().mark_explored(hdirection_);
      }
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

bool Algorithm::handle_abort_hdirection(WeightRef w)
{
  DoutEntering(dc::notice, "Algorithm::handle_abort_hdirection(" << w << ")");

  // Has a minimum been found at all?
  if (best_minimum_cubic_ == chain_.end())
  {
    Dout(dc::notice, "Aborting going " << hdirection_ << " and terminating search because no minimum has has been found at all!");
    return false;
  }

  // Jump back to the best minimum and continue in the opposite hdirection.
  Dout(dc::notice, "Aborting exploring " << hdirection_ << " of the minimum at " << best_minimum_cubic_->extreme_w() << ".");

  //
  // Go back in time to when we were at the best minimum so far, and then go into the opposite direction.
  //

  ASSERT(hdirection_ != HorizontalDirection::undecided);
  set_hdirection(opposite(hdirection_));        // Change hdirection.
  last_extreme_cubic_ = best_minimum_cubic_;
  saw_minimum_ = true;
  Dout(dc::notice, "Set last_extreme_cubic_ to [" << last_extreme_cubic_->label() << "].");
  cubic_used_ = chain_.end();

  return handle_local_extreme_jump(w);
}

bool Algorithm::handle_local_extreme_jump(WeightRef w)
{
  DoutEntering(dc::notice, "Algorithm::handle_local_extreme_jump(" << w << ")");

  // Set last_extreme_cubic_ to the local extreme that you want to jump (back) to.
  LocalExtreme const& local_extreme = last_extreme_cubic_->local_extreme();

  // Was this local extreme already explored in both directions?
  if (local_extreme.done())
  {
    // It shouldn't be possible that the best minimum wasn't explored into both directions yet, but an adjacent local extreme was.
    ASSERT(last_extreme_cubic_ == best_minimum_cubic_);
    handle_finish(w, {});
    return true;
  }

  // Was this local extreme already explored in hdirection?
  if (local_extreme.is_explored(hdirection_))
    return handle_abort_hdirection(w);

  Dout(dc::notice, "Jumping to " << local_extreme.label() << " with hdirection_ = " << hdirection_ << '.');

  // From here we want to go to an adjacent, hence opposite extreme.
  next_extreme_type_ = opposite(local_extreme.get_extreme_type());

  // Get neighbor in that direction, if any.
  SampleNode::const_iterator neighbor =
    (hdirection_ == HorizontalDirection::left) ? local_extreme.left_neighbor() : local_extreme.right_neighbor();

  // If this extreme has a neighbor then immediately jump to that instead.
  if (neighbor != chain_.end())
  {
    last_extreme_cubic_ = neighbor;
    Dout(dc::notice, "Set last_extreme_cubic_ to [" << last_extreme_cubic_->label() << "].");
    // We already saw a minimum.
    ASSERT(saw_minimum_);
    // Don't go back towards that minimum from this neighbor.
    neighbor->local_extreme().mark_explored(opposite(hdirection_));
    return handle_local_extreme_jump(w);
  }

  // Restore small_step_.
  small_step_ = last_extreme_cubic_->scale().value();

  if (local_extreme.opposite_direction_is_fourth_degree_extreme())
    handle_fourth_degree_approximation_jump(w, local_extreme.opposite_direction_w(), local_extreme.opposite_direction_Lw(), {});
  else
    handle_last_extreme_plus_step(w, {});

  // Now expected_Lw_ was updated. Don't forget to restore the energy too.
  energy_.set(best_minimum_energy_, expected_Lw_);

  // Determine the correct jump point (w might still not fall in the correct range).
  initialize_range(last_extreme_cubic_->extreme_w());
  move_into_range(w);

  return true;
}

bool Algorithm::handle_local_extreme(WeightRef w)
{
  DoutEntering(dc::notice, "Algorithm::handle_local_extreme(" << w << ")" << (state_ == IterationState::extra_sample ? " [Extra Sample]" : ""));
  RESTART

  auto extra_sample = state_ != IterationState::extra_sample ? chain_.end() : chain_.last();

  if (state_ == IterationState::extra_sample)
  {
    // The w value is still at that of the (requested) extra sample.
    ASSERT(w == extra_sample->w());

    // If we had only had two usable samples to fit a fourth degree polynomial,
    // then the third, extra sample should always have been added on the other
    // side of the sample close to the extreme (B) and thus never in the middle
    // of last_extreme_cubic_!
    //
    // That is, the situation should be one of:
    //
    // 1.   ------A*--------B-·-----------C---
    // 2.   ------A*----------·-B---------C---
    // 3.   ------C-----------·-B*--------A---
    // 4.   ------C---------B*·-----------A---
    //                        ^
    //                        |
    //                actual local extreme (next_extreme_type_)
    //
    // where C is the just added extra sample, and the node marked with a * is the last_extreme_cubic_.
    ASSERT(extra_sample == chain_.begin() || std::prev(extra_sample) != last_extreme_cubic_);
    // This assures that the last_extreme_cubic_ wasn't changed since the last entry.
  }
  else
  {
    // The 'local extreme cubic' may not contain an extreme of the opposite type: SampleNode has only one LocalExtreme pointer!
    // If this fails then we need to insert an extra sample to cut cubic_used_ in two.
    if (AI_UNLIKELY(cubic_used_->has_extreme(opposite(next_extreme_type_))))
    {
      // Apparently cubic_used_ contains two extremes. They should occur in the expected order.
      // If cubic_used_ == chain_.last() then we have case 4 (or 3), in which case the
      // next_extreme_type_ (B) must be on the left of the opposite extreme that is also part of B*-A.
      // Otherwise we have case 2 (or 1) and the next_extreme_type_ must be the extreme in A*-B that
      // is on the right.
      // If in fact we are dealing with case 1 or 3 then the current extreme type (next_extreme_type_)
      // is not contained in the cubic_used_ and has_extreme_order returns false in that case;
      // therefore begin with testing if the current extreme is also part of cubic_used_!
      ASSERT(!cubic_used_->has_extreme(next_extreme_type_) ||
          has_extreme_order(cubic_used_->type(), next_extreme_type_,
            cubic_used_ == chain_.last() ? HorizontalDirection::left : HorizontalDirection::right));
      // Request an extra sample in the middle of cubic_used_.
      handle_divide_last_extreme_cubic(w, {});
      return true;
    }

    // Remember the node containing the cubic that found this extreme.
    last_extreme_cubic_ = cubic_used_;
    Dout(dc::notice, "Set last_extreme_cubic_ to [" << last_extreme_cubic_->label() << "].");

    // Mark this sample as local extreme "cubic".
    last_extreme_cubic_->set_local_extreme(next_extreme_type_, chain_.end());
    LocalExtreme const& local_extreme = last_extreme_cubic_->local_extreme();

    // We should get here with a w value that was already set to the critical point of the latest cubic approximation.
    ASSERT(w == last_extreme_cubic_->extreme_w());

    // Because ascii-art is failing me big time, lets represent a local minimum like
    //
    //       \     /
    //        \   /
    //         `·´
    // as follows:
    //       ---·---
    //
    // where a · represents the local extreme and the rest of the function is simply flattened.
    //
    // A SampleNode (A or B, below) is marked as local extreme if its cubic has an extreme of
    // the desired type that is close to the last SampleNode (B, which is chain_.last()).
    //
    // We can distinguish the following cases (where last_extreme_cubic_ is marked with a * (the local extreme cubic)).
    //
    // 1.   ------A*--------B-·-----------C---          cubic: A-B has extreme in · close to B. A is marked as local extreme.
    //          <-^         ^->                         Start search to the left at A. Start search to the right at B.
    // 2.   ------A*----------·-B---------C---          cubic: A-B has extreme in · close to B. A is marked as local extreme.
    //          <-^             ^->                     Start search to the left at A. Start search to the right at B.
    // 3.   ------C-----------·-B*--------A---          cubic: B-A has extreme in · close to B. B is marked as local extreme.
    //          <-^             ^->                     Start search to the left at C. Start search to the right at B.
    // 4.   ------C---------B*·-----------A---          cubic: B-A has extreme in · close to B. B is marked as local extreme.
    //          <-^         ^->                         Start search to the left at C. Start search to the right at B.
    //
    // Note that in the first two cases, A is left of B and therefore contains the cubic.
    // Therefore A is last_extreme_cubic_, the one marked as local extreme (and B is chain_.last()).
    // In the last two cases, B is left of A and therefore B contains the cubic and is marked as local extreme:
    // B is both last_extreme_cubic_ and chain_.last().
    //
    // Let next_extreme_type_ be a ╿maximum│minimum╽, then
    //
    // in case 1. the type of A-B can have get_end(A-B, right) = ╿uphill│downhill╽
    //     and is not allowed to contain the opposite extreme therefore the type of A-B can be: ╿`/`│`\` or `^\`╽,
    //     in which case B-C is ╿`/\/`, `/\` or `/\_`│`\/`, `\/‾` or `\/\`╽ (get_end(B-C, left) is ╿uphill│downhill╽
    //     and B-C begins with a ╿maximum│minimum╽)
    //   or the type of A-B can have get_end(A-B, right) = ╿flat_high│flat_low╽
    //     and is not allowed to contain the opposite extreme therefore the type of A-B can only be ╿`/‾`│`\_`╽,
    //     in which case B-C is ╿`‾\`, `‾\_` or `‾\/`│`_/`, `_/‾` or `_/\`╽ (get_end(B-C, left) is ╿flat_high│flat_low╽ too
    //     and thus B-C begins with a ╿maximum│minimum╽).
    //
    // in case 2. the type of A-B can end on a ╿maximum│minimum╽
    //     and is not allowed to contain the opposite extreme therefore the type of A-B can only be ╿`/\`│`\/`╽,
    //     in which case get_end(B-C, left) is ╿downhill│uphill╽
    //     (the type begins with a ╿`\`: `\`, `\_`, `\/`, `\/‾` or `\/\`│`/`: `/`, `/^`, `/‾`, `/\/`, `/\` or `/\_`╽)
    //   or the type of A-B can have get_end(A-B, right) = ╿flat_high│flat_low╽
    //     and is not allowed to contain the opposite extreme therefore the type of A-B can only be ╿`/‾`│`\_`╽,
    //     in which case B-C is ╿`‾\`, `‾\_` or `‾\/`│`_/`, `_/‾` or `_/\`╽ (get_end(B-C, left) is ╿flat_high│flat_low╽ too
    //     and thus B-C begins with a ╿maximum│minimum╽).
    //
    // in case 3. the type of B-A can have get_end(B-A, left) = ╿downhill│uphill╽
    //     and is not allowed to contain the opposite extreme therefore the type of B-A can be: ╿`\`│`/` or `/^`╽,
    //     in which case C-B is ╿`/\`, `_/\` or `\/\`│`\/`, `‾\/` or `/\/`╽ (get_end(C-B, right) is ╿downhill│uphill╽
    //     and C-B ends on a ╿maximum│minimum╽)
    //   or the type of B-A can have get_end(B-A, left) = ╿flat_high│flat_low╽
    //     and is not allowed to contain the opposite extreme therefore the type of B-A can only be: ╿`‾\`│`_/`╽,
    //     in which case C-B is ╿`/‾`, `_/‾` or `\/‾`│`\_`, `‾\_` or `/\_`╽) (get_end(C-B, right) is ╿flat_high│flat_low╽ too
    //     and thus C-B ends on a ╿maximum│minimum╽).
    //
    // in case 4. the type of B-A can begin with a ╿maximum│minimum╽
    //     and is not allowed to contain the opposite extreme therefore the type of B-A can only be: ╿`/\`│`\/`╽,
    //     in which case get_end(C-B, right) is ╿uphill│downhill╽
    //     (the type ends on a ╿`/`: `/`, `_/`, `/\/`, `\/` or `‾\/`│`\`: `\`, `^\`, `‾\`, `\/\`, `/\` or `_/\`╽)
    //   or the type of B-A can have get_end(B-A, left) = ╿flat_high│flat_low╽
    //     and is not allowed to contain the opposite extreme therefore the type of B-A can only be: ╿`‾\`│`_/`╽,
    //     in which case C-B is ╿`/‾`, `_/‾` or `\/‾`│`\_`, `‾\_` or `/\_`╽) (get_end(C-B, right) is ╿flat_high│flat_low╽ too
    //     and thus C-B ends on a ╿maximum│minimum╽).
    //

    // Find in both directions the first cubic that has an extreme, and/or is marked as local extreme.

    // Case 3 or 4? Otherwise it's case 1 or 2.
    bool case34 = last_extreme_cubic_ == chain_.last();
    for (int d = -1; d <= 1; d += 2)
    {
      HorizontalDirection const direction = static_cast<HorizontalDirection>(d);

      // Determine the starting node (the node with the arrow in the ascii art).
      // 1.   ------A*--------B-·-----------C---
      //          <-^
      // 2.   ------A*----------·-B---------C---
      //          <-^
      // 3.   ------C-----------·-B*--------A---
      //                          ^->
      // 4.   ------C---------B*·-----------A---
      //                      ^->
      // First set start at the default (the node marked with a *).
      SampleNode::const_iterator start = last_extreme_cubic_;

      // However, if the direction is not left for case 1 and 2, or not right for case 3 and 4...
      if (case34 != (direction == HorizontalDirection::right))
      {
        if (!case34)
        {
          // Direction is right; start should have been B (if cubic B-C exists).
          // 1.   ------A*--------B-·-----------C---
          //                      ^->
          // 2.   ------A*----------·-B---------C---
          //                          ^->
          if (!start->is_cubic())       // Can't go further to the right?
            continue;                   // There is no opposite extreme on the right of A.
          ++start;
        }
        else
        {
          // Direction is left; start should have been C (if that exists).
          // 3.   ------C-----------·-B*--------A---
          //          <-^
          // 4.   ------C---------B*·-----------A---
          //          <-^
          if (start == chain_.begin())  // Can't go further to the left?
            continue;                   // There is no opposite extreme on the left of B.
          --start;
        }
      }

      SampleNode::const_iterator node = start;

      // Advance until we find a cubic X-Y with the a local extreme between X and Y,
      // or we run into another local extreme cubic (where X->is_local_extreme() is true,
      // which does not necessarily imply the former if the used extreme falls outside the range between X and Y).

      // The first time (while node == start) we only want to stop if node is a local extreme that is of the opposite type.
      if (node->is_cubic() &&
          (!node->is_local_extreme() || node->local_extreme().get_extreme_type() != opposite(next_extreme_type_)) &&
          !node->has_extreme(opposite(next_extreme_type_)))
      {
        if (direction == HorizontalDirection::left)
        {
          do
          {
            if (node == chain_.begin())
            {
              node = std::prev(chain_.end());   // We failed to find an opposite extreme; cause node->is_cubic() to be false.
              break;
            }
            --node;
          }
          while (!node->is_local_extreme() && !node->has_extreme());
        }
        else
        {
          do
          {
            ++node;
          }
          while (node->is_cubic() && !node->is_local_extreme() && !node->has_extreme());
        }
      }

#if CW_DEBUG
      if (node->is_cubic())
      {
        if (node == start)
        {
          // If node didn't change then it can't have been a local extreme cubic: if it was, then it can't also
          // be a local extreme cubic of the opposite type, nor is it allowed to contain an extreme of the opposite type.
          // If node didn't change and it not last_extreme_cubic_ then still it is not allowed have next_extreme_type_
          // as first local extreme, of course (extremes must alternate).
          ASSERT(start != last_extreme_cubic_);
        }
        else
        {
          // We found the first cubic that either has a local extreme, or is a local extreme cubic.
          // Since node was changed, it can impossibly contain the same extreme as where we came from,
          // unless it follows the opposite extreme in the same cubic.
          ASSERT(!node->has_extreme(next_extreme_type_) || has_extreme_order(node->type(), next_extreme_type_, direction));
        }

        // If the node does not have the sought for local extreme then it must be a local extreme cubic
        // that has its extreme outside of the cubic on the side of direction; aka we stopped one node too soon.
        // However, it is also possible that we skipped a local extreme, and we found a local extreme of
        // the same type as next_extreme_type_. In that case don't assert here, but assert below.
        ASSERT(node->has_extreme(opposite(next_extreme_type_)) ||
            (node->is_local_extreme() &&
             ((direction == HorizontalDirection::left && node == chain_.begin()) ||
              (direction == HorizontalDirection::right && !std::next(node)->is_cubic()) ||
              (direction == HorizontalDirection::left ? std::prev(node) : std::next(node))->has_extreme(opposite(next_extreme_type_)) ||
              // Do not assert here, but assert below.
              node->local_extreme().get_extreme_type() != opposite(next_extreme_type_))));
      }
#endif

      if (node->is_cubic() && !node->is_local_extreme())
      {
        // If this is not a local extreme, then the next one must be - or this local extreme wasn't found yet.
        if (direction == HorizontalDirection::left)
        {
          if (node != chain_.begin())
            --node;
        }
        else
          ++node;
      }

      if (node->is_cubic() && node->is_local_extreme())
      {
        // If we found a neighboring local extreme, it better be of the opposite type.
        ASSERT(node->local_extreme().get_extreme_type() == opposite(next_extreme_type_));
        if (direction == HorizontalDirection::left)
        {
          last_extreme_cubic_->local_extreme().set_left_neighbor(node);
          node->local_extreme().set_right_neighbor(last_extreme_cubic_);
        }
        else
        {
          last_extreme_cubic_->local_extreme().set_right_neighbor(node);
          node->local_extreme().set_left_neighbor(last_extreme_cubic_);
        }
      }

    } // Another direction?

#ifdef CWDEBUG
    // Draw the local extreme sample.
    event_server_.trigger(AlgorithmEventType{new_local_extreme_event, *last_extreme_cubic_, std::to_string(last_extreme_cubic_->label())});
    event_server_.trigger(AlgorithmEventType{scale_erase_event});
#endif

    // Get neighbor in the direction that we came from, if any.
    SampleNode::const_iterator neighbor =
      (hdirection_ == HorizontalDirection::right) ? local_extreme.left_neighbor() : local_extreme.right_neighbor();

    if (neighbor != chain_.end())
    {
      // If we find a new local extreme by exploring in the direction left/right from a previous extreme and it is adjacent
      // (there is no extreme in between), then the previous extreme is marked as Explored in that direction.
      neighbor->local_extreme().mark_explored(hdirection_);
      //FIXME: uncomment his and remove the ASSERT below.
      //saw_minimum_ = false;
    }

    // If we already saw a minimum, going from adjacent local extremes to the next, then also mark the newly found extreme
    // as Explored in the opposite direction (towards that already seen extreme).
    if (saw_minimum_)
    {
      // If saw_minimum_ is true, then we must have a neighbor.
      // If this fails we must have skipped local extremes between the one that caused saw_minimum_ to be set and this one.
      // FIXME: in the end this should be allowed. Uncomment the line above.
      ASSERT(neighbor != chain_.end());
      local_extreme.mark_explored(opposite(hdirection_));

      // Did we just find a minimum?
      if (next_extreme_type_ == ExtremeType::minimum)
      {
        // Find the neighbor of the neighbor, which has to be the minimum that we saw (aka, it must exist).
        SampleNode::const_iterator neighbor_of_neighbor =
          (hdirection_ == HorizontalDirection::right) ? neighbor->local_extreme().left_neighbor() : neighbor->local_extreme().right_neighbor();
        ASSERT(neighbor_of_neighbor != chain_.end());
        neighbor_of_neighbor->local_extreme().mark_explored(hdirection_);
      }
    }

    // Keep track of the best minimum so far; or abort if this minimum isn't better than one found before.
    if (next_extreme_type_ == ExtremeType::minimum)
    {
      Dout(dc::notice, "new minimum = {" << last_extreme_cubic_->extreme_w() << ", " <<
          last_extreme_cubic_->local_extreme().extreme_Lw() << "}");
      // We were looking for a minimum.
      ASSERT(last_extreme_cubic_->local_extreme().get_extreme_type() == ExtremeType::minimum);
      if (best_minimum_cubic_ == chain_.end() ||
          best_minimum_cubic_->local_extreme().extreme_Lw() > last_extreme_cubic_->local_extreme().extreme_Lw())
      {
        best_minimum_cubic_ = last_extreme_cubic_;
        Dout(dc::notice, "best_minimum_cubic_ set to " << *best_minimum_cubic_);
        best_minimum_energy_ = energy_.energy();
      }
      if (last_extreme_cubic_ != best_minimum_cubic_)
      {
        // The new minimum isn't better than what we found already. Stop going into this direction.
        Dout(dc::notice, "The new minimum isn't better than what we found already. "
            "Stop going into the direction " << hdirection_ << ".");
        return false;
      }
    }

    // Set saw_minimum_ if this is a minimum.
    if (next_extreme_type_ == ExtremeType::minimum)
      saw_minimum_ = true;
  }

  double const extreme_w = last_extreme_cubic_->extreme_w();
  double const scale_value = last_extreme_cubic_->scale().value();

  // Update small_step_ (value() returns a positive value).
  small_step_ = scale_value;
  Dout(dc::notice, "small_step_ set to " << small_step_);

  // Find all samples around last_extreme_ that are within scale.
  SampleNode::const_iterator left_edge = last_extreme_cubic_;
  SampleNode::const_iterator right_edge = left_edge->next_node();

  constexpr int max_samples_per_side = 3;
  std::array<SampleNode::const_iterator, 2 * max_samples_per_side + 1> usable_samples;

  int number_of_usable_samples = 0;
  int local_extreme_index = 0;          // Index into usable_samples for the SampleNode that is closest to extreme_w.
  double smallest_dist = std::numeric_limits<double>::max();
  auto add_sample = [&, extreme_w](SampleNode::const_iterator node){
    if (std::abs(node->w() - extreme_w) < smallest_dist)
    {
      smallest_dist = std::abs(node->w() - extreme_w);
      local_extreme_index = number_of_usable_samples;
    }
    usable_samples[number_of_usable_samples++] = node;
  };
  // Always add left_edge and right_edge.
  add_sample(left_edge);
  add_sample(right_edge);
  bool left_edge_is_on_the_left_of_extreme = left_edge->w() < extreme_w;
  bool right_edge_is_on_the_left_of_extreme = right_edge->w() < extreme_w;
  int remaining_left_samples = max_samples_per_side - left_edge_is_on_the_left_of_extreme - right_edge_is_on_the_left_of_extreme;
  int remaining_right_samples = max_samples_per_side - !left_edge_is_on_the_left_of_extreme - !right_edge_is_on_the_left_of_extreme;

  double const left_scale = last_extreme_cubic_->scale().left_edge_w();       // The smallest w value for which the cubic is still matching according to scale.
  bool saw_extra_sample = false;
  while (left_edge != chain_.begin() && remaining_left_samples > 0)
  {
    double const prev_w = left_edge->w();
    --left_edge;
    // Always add extra_sample if it exists.
    bool is_extra_sample = left_edge == extra_sample;
    saw_extra_sample |= is_extra_sample;
    if (!is_extra_sample)
    {
      // Stop if we went too far to the left.
      if (left_edge->w() < left_scale)
        break;
      // Calculate the distance to the previous sample, if any.
      double dist = prev_w - left_edge->w();
      ASSERT(dist > 0.0);
      if (dist < 0.001 * scale_value)
        continue;
    }
    add_sample(left_edge);
    --remaining_left_samples;
  }
  double const right_scale = last_extreme_cubic_->scale().right_edge_w();       // The largest w value for which the cubic is still matching according to scale.
  double prev_w = right_edge->w();
  while (++right_edge != chain_.end() && remaining_right_samples > 0)
  {
    // Always add extra_sample if it exists.
    bool is_extra_sample = right_edge == extra_sample;
    saw_extra_sample |= is_extra_sample;
    if (!is_extra_sample)
    {
      // Stop if we went too far to the right.
      if (right_edge->w() > right_scale)
        break;
      // Calculate the distance to the previous sample, if any.
      double dist = right_edge->w() - prev_w;
      ASSERT(dist > 0.0);
      if (dist < 0.001 * scale_value)
        continue;
    }
    add_sample(right_edge);
    --remaining_right_samples;
    prev_w = right_edge->w();
  }
  // Always add extra_sample.
  if (!saw_extra_sample && extra_sample != chain_.end())
  {
    // We added at most twice max_samples_per_side.
    ASSERT(number_of_usable_samples < usable_samples.size());
    add_sample(extra_sample);
  }

#ifdef CWDEBUG
  {
    Dout(dc::notice, "Number of samples within scale range: " << number_of_usable_samples);
    auto sorted_usable_samples = usable_samples;
    std::sort(sorted_usable_samples.begin(), sorted_usable_samples.begin() + number_of_usable_samples,
        [](SampleNode::const_iterator node1, SampleNode::const_iterator node2){ return node1->w() < node2->w(); });
    for (int i = 0; i < number_of_usable_samples; ++i)
      Dout(dc::notice, *sorted_usable_samples[i]);
  }
#endif

  // There must be at least one usable sample: one of the samples used to fit last_extreme_cubic_ may
  // lay further away than 1.1 times the scale, but not both.
  ASSERT(number_of_usable_samples > 0);

  // After asking for one extra sample we should now have two or three samples.
  ASSERT(state_ != IterationState::extra_sample || number_of_usable_samples >= 2);

  // If we do not already have at least three relevant samples then get another one.
  if (number_of_usable_samples < 3)
  {
    ASSERT(number_of_usable_samples == 2);

    // Find the sample with the greatest distant to the local extreme.
    double max_offset = 0.0;
    for (int i = 0; i < number_of_usable_samples; ++i)
    {
      double offset = usable_samples[i]->w() - extreme_w;
      if (std::abs(offset) > std::abs(max_offset))
        max_offset = offset;
    }

    // Request an extra sample on the other side of the critical point than that we already have.
    HorizontalDirection direction = max_offset > 0.0 ? HorizontalDirection::left : HorizontalDirection::right;
    handle_get_extra_sample(w, extreme_w + last_extreme_cubic_->scale().step(direction), {});

    // Remember this direction in hdirection_. This is used when we can't find any local minimum in the fourth degree polynomial.
    if (hdirection_ == HorizontalDirection::undecided)
      set_hdirection(direction);

    return true;
  }

  // After finding a maximum we want to find a minimum and visa versa. Change next_extreme_type_.
  next_extreme_type_ = opposite(next_extreme_type_);
  Dout(dc::notice, "next_extreme_type_ is toggled to " << next_extreme_type_ << ".");

  // Use the w value of the Sample that is closest to the found extreme.
  Sample const& w2 = *usable_samples[local_extreme_index];
  double const w2_1 = w2.w();

  // Brute force find the two samples that, together with the local extreme sample w2, have the largest spread.
  int i0 = local_extreme_index == 0 ? 1 : 0;
  int i1 = local_extreme_index == i0 + 1 ? i0 + 2 : i0 + 1;
  double best_spread = 0.0;
  if (number_of_usable_samples > 3)
  {
    for (int t0 = i0; t0 < number_of_usable_samples; ++t0)
    {
      if (t0 == local_extreme_index)
        continue;
      for (int t1 = t0 + 1; t1 < number_of_usable_samples; ++t1)
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

  math::Polynomial derivative = fourth_degree_approximation.derivative();
  Dout(dc::notice, "derivative = " << derivative);

#ifdef CWDEBUG
//  event_server_.trigger(AlgorithmEventType{derivative_event, derivative});
#endif

  double remainder;
  auto quotient = derivative.long_division(extreme_w, remainder);
  Dout(dc::notice, "quotient = " << quotient << " with remainder " << remainder);

  auto second_derivative = quotient.derivative();
  Dout(dc::notice, "second_derivative = " << second_derivative);

  // Remember the left- and right-most samples that were used to fit the fourth degree polynomial.
  last_extreme_cubic_->local_extreme().set_edge_nodes(fourth_degree_approximation, usable_samples, local_extreme_index, i0, i1);

#ifdef CWDEBUG
//  event_server_.trigger(AlgorithmEventType{quotient_event, quotient});
#endif

  std::array<double, 2> quotient_roots;
  int number_of_quotient_roots = quotient.get_roots(quotient_roots);
  ASSERT(number_of_quotient_roots == 0 || number_of_quotient_roots == 2);

  // The following five possibilities exist:
  //
  // 1) There are no roots:
  //        We need to do a step in the right direction, and set opposite_direction_w to a step in the opposite direction.
  // 2) There are two roots, one of which is of the wrong type.
  //    This can only happen if the found extreme is on the edge, for example: <minimum> <maximum> <irrelevant minimum>:
  //    Then, either                                                     extreme_w-^         ^-the two roots-^
  //    2a) the relevant root is on the wrong side (hdirection_ is already set to left or right):
  //        We need to do a step in the direction of hdirection_, and set opposite_direction_w to the "relevant" root.
  //    2b) the relevant root is on the correct side (possibly because hdirection_ is still undecided):
  //        We need to jump to the relevant root, and set opposite_direction_w to a step in the opposite direction.
  // 3) There are two roots that both are usable (hdirection_ is undecided).
  //        We need to jump to the lowest root, and set opposite_direction_w to the other root.
  // 4) There are two roots, one of which is on the wrong side (hdirection_ is already set to left or right).
  //        We need to jump to the root that is on the correct side, and set opposite_direction_w to the other root.

  bool opposite_direction_is_fourth_degree_extreme;
  double opposite_direction_w, opposite_direction_Lw;
  if (number_of_quotient_roots == 0)
  {
    // Case 1.
    //   We need to do a step in the right direction, and set opposite_direction_w to a step in the opposite direction.

    if (hdirection_ == HorizontalDirection::undecided)
      set_hdirection(last_extreme_cubic_->scale().direction_trend());
    handle_last_extreme_plus_step(w, {});

    // Must call handle_last_extreme_plus_step(w, {}) for the opposite direction.
    opposite_direction_is_fourth_degree_extreme = false;

    Dout(dc::notice, "with no other usable local extrema! Using step in both directions.");

    // Continue finding the next extreme.
    state_ = IterationState::find_extreme;
  }
  else if (extreme_w > quotient_roots[0] == extreme_w > quotient_roots[1])      // Both roots are on the same side as the found extreme?
  {
    // Case 2.
    // The root closest to extreme_w has the correct type (opposite of the one we just found).
    int relevant_root = std::abs(extreme_w - quotient_roots[0]) < std::abs(extreme_w - quotient_roots[1]) ? 0 : 1;
    if (hdirection_ != HorizontalDirection::undecided && wrong_direction(quotient_roots[relevant_root] - extreme_w))
    {
      // Case 2a.
      //   We need to do a step in the direction of hdirection_, and set opposite_direction_w to the "relevant" root.
      handle_last_extreme_plus_step(w, {});

      // Must call handle_fourth_degree_approximation_jump(w, quotient_roots[relevant_root],
      // fourth_degree_approximation(opposite_direction_w), {}) for the opposite direction.
      opposite_direction_is_fourth_degree_extreme = true;
      opposite_direction_w = quotient_roots[relevant_root];
      opposite_direction_Lw = fourth_degree_approximation(opposite_direction_w);

      Dout(dc::notice, "with one other usable local extreme, but on the wrong side; using a step into " << hdirection_ <<
          " -- while using that root for opposite_direction_w: " << opposite_direction_w);

      // Continue finding the next extreme.
      state_ = IterationState::find_extreme;
    }
    else
    {
      // Case 2b.
      //   We need to jump to the relevant root, and set opposite_direction_w to a step in the opposite direction.
      handle_fourth_degree_approximation_jump(w, quotient_roots[relevant_root], fourth_degree_approximation(quotient_roots[relevant_root]), {});

      if (hdirection_ == HorizontalDirection::undecided)
        set_hdirection(w < w2_1 ? HorizontalDirection::left : HorizontalDirection::right);

      auto& scale = last_extreme_cubic_->scale();
      // Must call handle_last_extreme_plus_step(w, {}) for the opposite direction.
      opposite_direction_is_fourth_degree_extreme = false;

      Dout(dc::notice, "with one other usable local extreme at " << w << ", which is in the correct direction (" << hdirection_ << ")" <<
          " -- we will do a scale step in the opposite direction.");

      state_ = IterationState::backtracking;
    }
  }
  else
  {
    std::array<double, 2> expected_Lw;
    for (int quotient_root_index = 0; quotient_root_index < 2; ++quotient_root_index)
      expected_Lw[quotient_root_index] = fourth_degree_approximation(quotient_roots[quotient_root_index]);
    int best_root;
    if (hdirection_ == HorizontalDirection::undecided)
    {
      // Case 3.
      //   We need to jump to the lowest root, and set opposite_direction_w to the other root.
      best_root = expected_Lw[0] < expected_Lw[1] ? 0 : 1;

      Dout(dc::notice|continued_cf, "with two other usable local extrema; the lowest one being at " << quotient_roots[best_root] << " -- ");
    }
    else
    {
      // Case 4.
      //   We need to jump to the root that is on the correct side, and set opposite_direction_w to the other root.
      best_root = wrong_direction(quotient_roots[1] - extreme_w) ? 0 : 1;

      Dout(dc::notice|continued_cf, "with two other usable local extrema; one in hdirection (" << hdirection_ << ") at " <<
          quotient_roots[best_root] << " -- ");
    }
    handle_fourth_degree_approximation_jump(w, quotient_roots[best_root], expected_Lw[best_root], {});

    // Must call handle_fourth_degree_approximation_jump(w, quotient_roots[1 - best_root], expected_Lw[1 - best_root], {}) for
    // the opposite direction.
    opposite_direction_is_fourth_degree_extreme = true;
    opposite_direction_w = quotient_roots[1 - best_root];
    opposite_direction_Lw = expected_Lw[1 - best_root];

    Dout(dc::finish, "and using the other one for opposite_direction_w: " << opposite_direction_w);

    state_ = IterationState::backtracking;
  }

  // Also remember the opposite direction (in case we jump back here).
  last_extreme_cubic_->local_extreme().set_opposite_direction(
      opposite_direction_is_fourth_degree_extreme, opposite_direction_w, opposite_direction_Lw);

  if (hdirection_ == HorizontalDirection::undecided)
  {
    // Now that the decision on which hdirection_ we explore is taken, store that decision.
    set_hdirection(w < w2_1 ? HorizontalDirection::left : HorizontalDirection::right);
  }

  initialize_range(extreme_w);
  move_into_range(w);

  return true;
}

void Algorithm::initialize_range(double extreme_w)
{
  DoutEntering(dc::notice, "Algorithm::initialize_range(" << extreme_w << ")");

  // hdirection_ must be set before calling this function.
  ASSERT(hdirection_ != HorizontalDirection::undecided);

  using enum CubicToNextSampleType;

  // Determine the correct range, by setting left_of_ and right_of_, given that
  // we will be looking for next_extreme_type_ in the direction hdirection_.
  left_of_ = chain_.end();
  right_of_ = chain_.end();

  // Get the node immediately on the left and right of extreme_w.
  chain_.find_larger(extreme_w);
  SampleNode::const_iterator right_node = chain_.larger();
  SampleNode::const_iterator left_node = right_node == chain_.begin() ? chain_.end() : std::prev(right_node);

#ifdef CWDEBUG
  // Note that extreme_w is an extreme of the opposite type of next_extreme_type_.
  if (left_node != chain_.end() && left_node->is_cubic())    // Can't check the cubic type if it doesn't exist.
  {
    if (next_extreme_type_ == ExtremeType::maximum)
    {
      // This means that extreme_w is a minimum.
      ASSERT(has_minimum(left_node->type()));
      // Aka, the type is _/, \_, ‾\_, _/‾, \/, ‾\/, \/‾, \/\, /\/, /\_ or _/\.

      // If the first cubic also contains a maximum, then it must be of the type
      // ‾\_, _/‾, ‾\/, \/‾, \/\, /\/, /\_ or _/\.
    }
    else
    {
      // This means that extreme_w is a maximum.
      ASSERT(has_maximum(left_node->type()));
      // Aka, the type is ‾\, /‾, _/‾, ‾\_, /\, _/\, /\_, /\/, \/\, \/‾ or ‾\/.

      // If the first cubic also contains a minimum, then it must be of the type
      // _/‾, ‾\_, _/\, /\_, /\/, \/\ or ‾\/.
    }
  }
#endif

  if (hdirection_ == HorizontalDirection::right)
  {
    // For example, when looking for a maximum (next_extreme_type_) and the situation is:
    //
    // last_extreme_cubic_
    //     |           .--- the maximum that we want to find next (it's on the right);
    //     |           |    it is on the right of left_node, NOT on the right of right_node.
    //     |           v
    //    \v          .-.
    //     O         /   \
    //      \       /     O
    //       \     /      ^\
    //        O   /       |
    //         `·´   right_node (the first node larger than extreme_w)
    //        ^ ^
    //        | |
    //        | extreme_w = last_extreme_cubic_->extreme_w()
    //        |
    //   left_node (the node immediately on the left of right_node)
    //
    // We have left_node point to a cubic with two extremes, and need right_of_ to be set to left_node (and not advance it),
    // because the maximum that we're looking for IS the maximum between left_node and right_node.
    //
    // Note that it is possible that right_node doesn't exist:
    //
    // last_extreme_cubic_
    //     |           .--- the maximum that we want to find next (it's on the right);
    //     |           |    it is on the right of left_node.
    //     |           v
    //    \v          .-.
    //     O         /   \
    //      \       /     \
    //       \     /       \
    //        O   /
    //         `·´   right_node = chain_.end() (there is no sample on the right of extreme_w).
    //        ^ ^
    //        | |
    //        | extreme_w = last_extreme_cubic_->extreme_w()
    //        |
    //   left_node (still std::prev(right_node): it's the last node in the chain)
    //
    // In which case we want right_of_ to be set to left_node, so that we'll look for a possible maximum beyond the last sample.
    //
    // The situation can also be the following:
    //
    //              last_extreme_cubic_
    //                   |     .--- the maximum that we want to find next (it's on the right of extreme_w);
    //                   |     |    it is on the right of right_node!
    //                   |     v
    //          M        |    .-.
    //         .-.       |   /   \
    //        /   \      |  O
    //       O     \     | /^
    //      /^      \    v/ |
    //     / |       \   O  afterwards we'll advance right_of_ to this sample (for example).
    //       |        `·´^
    //       |         ^ |
    //       |         | right_node (the first node larger than extreme_w)
    //       |         |
    //       |         extreme_w = last_extreme_cubic_->extreme_w()
    //       |
    //   left_node (the node immediately on the left of right_node)
    //
    // We have again left_node point to a cubic with two extremes, but need right_of_ to be set to right_node
    // because it can't be left_node: that would cause us to find the maximum on the left of extreme_w (M).
    //
    // Finally, left_node might not exist:
    //
    //              last_extreme_cubic_
    //                   |     .--- the maximum that we want to find next (it's on the right of extreme_w);
    //                   |     |    it is on the right of right_node!
    //                   |     v
    //          M        |    .-.
    //         .-.       |   /   \
    //        /   \      |  O
    //             \     | /^
    //              \    v/ |
    //               \   O  afterwards we'll advance right_of_ to this sample (for example).
    //                `·´^
    //                 ^ |
    //                 | right_node (the first node larger than extreme_w)
    //                 |
    //                 extreme_w = last_extreme_cubic_->extreme_w()
    //
    //   left_node = chain_.end() (there is no sample on the left of extreme_w).
    //
    // Also here we need right_of_ to be set to right_node.
    //
    right_of_ =
      (left_node != chain_.end() &&
       (!left_node->is_cubic() || has_extreme_order(left_node->type(), next_extreme_type_, HorizontalDirection::right))) ? left_node
                                                                                                                         : right_node;
    // Now that is settled, advance right_of_ as long as possible.
    while (right_of_->is_cubic() && !right_of_->has_extreme(next_extreme_type_))
      ++right_of_;

    if (right_of_->is_cubic())
      left_of_ = right_of_->next_node();
  }
  else
  {
    // Likewise, still assuming we're looking for a maximum (next_extreme_type_) and the situation is:
    //
    //       .---- the maximum that we want to find next (it's on the right of extreme_w);
    //       |     it is on the left of right_node, NOT on the left of left_node.
    //       v
    //      .-.           .-.
    //     /   \         O
    //    O     \       /
    //   /^      \     /
    //  / |       \   O <-- last_extreme_cubic_
    //    |        `·´^
    //    |         ^ `-- right_node (the first node larger than extreme_w)
    //    |         |
    //    |         extreme_w = last_extreme_cubic_->extreme_w()
    // left_node (the node immediately on the left of right_node)
    //
    // We have left_node point to a cubic with two extremes, and need left_of_ to be set to right_node (and not advance it),
    // because the maximum that we're looking for IS the maximum between left_node and right_node.
    //
    // Note that it is possible that left_node doesn't exist:
    //
    //       .---- the maximum that we want to find next (it's on the left of extreme_w);
    //       |     it is on the left of right_node.
    //       v
    //      .-.           .-.
    //     /   \         O
    //    /     \       /
    //   /       \     /
    //            \   O <-- last_extreme_cubic_
    //             `·´^
    //              ^ `-- right_node (the first node larger than extreme_w)
    //              |
    //              extreme_w = last_extreme_cubic_->extreme_w()
    // left_node = chain_.end() (there is no sample on the left of extreme_w).
    //
    // In which case we want left_of_ to be set to right_node, so that we'll look for a possible maximum beyond the last sample.
    //
    // The situation can also be the following:
    //
    // last_extreme_cubic_
    //  .--|--- the maximum that we want to find next (it's on the left of extreme_w);
    //  v  |    it is on the left of left_node.
    // .-. |           M
    //    \v          .-.
    //     O         /   \
    //      \       /     O
    //       \     /      ^\
    //        O   /       |
    //         `·´   right_node (the first node larger than extreme_w)
    //        ^ ^
    //        | |
    //        | extreme_w = last_extreme_cubic_->extreme_w()
    //        |
    //   left_node (the node immediately on the left of right_node)
    //
    // We have left_node point to a cubic with two extremes, but need left_of_ to be set to left_node
    // because it can't be right_node: that would cause us to find the maximum on the right of extreme_w (M).
    //
    // Finally, right_node might not exist:
    //
    // last_extreme_cubic_
    //  .--|--- the maximum that we want to find next (it's on the left of extreme_w);
    //  v  |    it is on the left of left_node.
    // .-. |           M
    //    \v          .-.
    //     O         /   \
    //      \       /
    //       \     /
    //        O   /
    //         `·´   right_node = chain_.end() (there is no sample on the right of extreme_w).
    //        ^ ^
    //        | |
    //        | extreme_w = last_extreme_cubic_->extreme_w()
    //        |
    //   left_node (the node immediately on the left of right_node)
    //
    // Also here we need left_of_ to be set to left_node.
    //
    // Note that in this case we want to use right_node if left_node doesn't exist (left_node == chain_.end()),
    // but left_node if right_node doesn't exist (left_node->is_cubic() is false).
    left_of_ =
      (left_node != chain_.end() &&
       (!left_node->is_cubic() || has_extreme_order(left_node->type(), next_extreme_type_, HorizontalDirection::right))) ? left_node
                                                                                                                         : right_node;
    while (left_of_ != chain_.begin() && !std::prev(left_of_)->has_extreme(next_extreme_type_))
      --left_of_;

    if (left_of_ != chain_.begin())
      right_of_ = std::prev(left_of_);
  }

  Dout(dc::notice, "New range: " << utils::print_using(*this, &Algorithm::print_range_on));

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{left_of_right_of_event,
      left_of_ == chain_.end() ? nullptr : &*left_of_,
      right_of_ == chain_.end() ? nullptr : &*right_of_});
#endif
}

#ifdef CWDEBUG
void Algorithm::print_range_on(std::ostream& os) const
{
  if (right_of_ == chain_.end())
    os << "(-inf";
  else
    os << "[[" << right_of_->label() << "]";
  if (left_of_ == chain_.end())
    os << ", +inf)";
  else
  {
    if (right_of_ != chain_.end())
      os << ' ' << to_utf8_art(left_of_->type()) << ' ';
    else
      os << ", ";
    os << "[" << left_of_->label() << "]]";
  }
}
#endif

} // namespace gradient_descent
