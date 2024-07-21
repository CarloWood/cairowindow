#include "sys.h"
#include "Algorithm.h"
#include <Eigen/Dense>
#include <numeric>
#ifdef CWDEBUG
#include "utils/print_using.h"
#endif

namespace gradient_descent {

bool Algorithm::operator()(Weight& w, double Lw, double dLdw)
{
  DoutEntering(dc::notice, "Algorithm::operator()(" << w << ", " << Lw << ", " << dLdw << ")");
  Dout(dc::notice, "hdirection_ = " << hdirection_ << "; next_extreme_type_ = " << next_extreme_type_ << "; state_ = " << state_);

#ifdef CWDEBUG
  // Erase all previous curves (if they exist).
  event_server_.trigger(AlgorithmEventType{reset_event});
  event_server_.trigger(AlgorithmEventType{difference_event, w, expected_Lw_, Lw});
#endif

  // If the new sample (w) is too close to the previous sample (Scale::negligible returns true)
  // then the new sample replaces the previous sample, current_is_replacement is set to true
  // and no new sample is added to the history.
  //
  // Otherwise, the new sample is added to the history and current_is_replacement is set to false.
  bool current_is_replacement;
  history_.add(w, Lw, dLdw, approximation_ptr_->scale(), current_is_replacement);

  // This function should never return a value whose difference with the previous sample is negligible
  // if there is only a single relevant sample in the history.
  ASSERT(!current_is_replacement || history_.relevant_samples() > 1);

  bool have_zero = false;
  Weight w_zero(w);
  if (state_ == IterationState::need_extra_sample)
  {
    // We're handling a local extreme, even though this is an extra sample.
    ASSERT(approximation_ptr_->is_extreme());
    // We already found a local extreme and history_.current() is that extreme.
    ASSERT(history_.total_number_of_samples() > 0);
    // w should not be negligible (handle_local_extreme was responsible for that).
    ASSERT(!current_is_replacement);

    // Fit a fourth degree polynomial through the local extreme and return the next point to jump to based on that.
    // This call sets w_zero and returns true, or returns false if the fitted fourth degree polynomial does not have usable extremes.
    have_zero = handle_local_extreme(w_zero);

    // Fall-through to update_approximation; which we want to call update_local_extreme_scale.
    ASSERT(approximation_ptr_ != &current_approximation_);
  }

  // While back tracking we're not really moving forward, we're just probing.
  // Therefore we shouldn't decrease the kinetic energy as if we traveled a distance.
  if (state_ != IterationState::back_tracking)
  {
    // Update kinetic energy. Returns false if too much energy was used.
    if (!update_energy())
      return handle_abort_hdirection(w);
  }

  // Handle the case where this sample is a local extreme.
  if (state_ == IterationState::local_extreme)
  {
    // This state is set when the curve that we're trying to find the extreme of locally looks like a cubic,
    // and the last sample is in the extreme that we're looking for (i.e. has a derivative close to zero).
    // Returns false if this exreme is a minimum but isn't better than the previously found best minimum.
    // Otherwise returns true and sets w to the new jump point.
    if (!handle_local_extreme(w))
      return handle_abort_hdirection(w);
    Dout(dc::notice, "Returning: " << std::setprecision(std::numeric_limits<double>::max_digits10) << w);
    return true;        // w was successfully updated by handle_local_extreme.
  }

  // last_region_ will be changed by the call to update_approximation, which (might) call find_extreme.
  bool const first_call = last_region_ == Region::unknown;
  // If first_call is true, then approximation_ptr_ shouldn't point to a LocalExtreme.
  ASSERT(!first_call || approximation_ptr_ == &current_approximation_);

  // Create/update a cubic approximation from this and the previous sample (or a line if this is the first sample).
  // This might return nonsense. new_w is only set if first_call is true, the Approximation is not part
  // of a LocalExtreme and the approximation has two relevant samples.
  double new_w = update_approximation(current_is_replacement);

  if (approximation_ptr_->number_of_relevant_samples() == 1)
  {
    // This function decides where to probe for a new sample next, based on a single sample, and changes w accordingly.
    handle_single_sample(w);
    Dout(dc::notice, "Returning: " << std::setprecision(std::numeric_limits<double>::max_digits10) << w);
    return true;
  }
  // Here we have two samples, therefore now new_w must be valid if first_call is true.

  // If have_zero is set then the current sample is an extra sample that was needed to
  // fit a fourth degree polynomial through a found local extreme, and that fit had
  // a usable extreme to jump to; hence jump there now. However, if the jump by
  // coincidence ends up very close to where we are already anyway, then there is no
  // need to return and ask for a new sample. In that case we might as well pretend
  // we already jumped and continue with the current sample.
  if (have_zero && !approximation_ptr_->scale().negligible(w - w_zero))
  {
    w = w_zero;
    state_ = IterationState::done;
    Dout(dc::notice, "Returning: " << std::setprecision(std::numeric_limits<double>::max_digits10) << w);
    return true;
  }

  // Decide where to jump next based on the current approximation (this call changes w).
  handle_approximation(w, first_call, new_w);

  if (state_ == IterationState::local_extreme)
  {
    // Do not return a value that would replace the current sample.
    if (approximation_ptr_->scale().negligible(w - history_.current().w()))
    {
      // Instead handle the local extreme immediately (this call changes w).
      if (!handle_local_extreme(w))
        return handle_abort_hdirection(w);
    }
  }

  Dout(dc::notice, "Returning: " << std::setprecision(std::numeric_limits<double>::max_digits10) << w);
  return true;
}

bool Algorithm::update_energy()
{
  // Get the height reached.
  double Lw = history_.current().Lw();

  if (state_ == IterationState::check_energy)
  {
    // Update the current kinetic energy. If this is an overshoot, abort this horizontal direction.
    //
    // This state is set when locally the curve looks like a cubic with its critical points
    // on the side that we just came from: in this case we move away from the all critical points
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

void Algorithm::reset_history()
{
  DoutEntering(dc::notice, "Algorithm::reset_history()");

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{scale_erase_event});
#endif
  // Use the Approximation object on the stack again (as opposed to one from a LocalExtreme).
  approximation_ptr_ = &current_approximation_;
  // Reset the cubic approximation.
  approximation_ptr_->reset();
  history_.reset();
}

#ifdef CWDEBUG
void Algorithm::set_algorithm_str(double new_w, char const* algorithm_str)
{
  algorithm_str_ = algorithm_str;
  Dout(dc::notice, std::setprecision(12) << history_.current().w() << " --> " << history_.total_number_of_samples() << ": " <<
      new_w << " [" << algorithm_str_ << "] [expected_Lw: " << expected_Lw_ << "]");
}
#endif

bool Algorithm::handle_local_extreme(Weight& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_local_extreme(" << w << ")");

  // Following adding the first extreme, we should have decided on a horizontal direction.
  ASSERT(hdirection_ != HorizontalDirection::undecided || extremes_.empty());

  if (state_ == IterationState::local_extreme)
  {
    auto prev_extreme = last_extreme_;

    // This means that history_.current() is a local extreme. Store it as an extreme.
    last_extreme_ =
      extremes_.emplace(hdirection_ == HorizontalDirection::right ? extremes_.end() : extremes_.begin(),
          history_.current(), std::move(*approximation_ptr_), energy_.energy());

    if (prev_extreme != extremes_.end())
    {
      // Add pointer back to the previous local extreme.
      last_extreme_->set_neighbor(opposite(hdirection_), &*prev_extreme);
    }

#ifdef CWDEBUG
    event_server_.trigger(AlgorithmEventType{new_local_extreme_event, last_extreme_->cp_sample(),
        std::to_string(history_.total_number_of_samples() - 1)});
    event_server_.trigger(AlgorithmEventType{scale_erase_event});
#endif
    // Switch approximation_ptr to the approximation stored in this extreme:
    // we need to keep updating it when new samples are added that match the same cubic.
    approximation_ptr_ = &last_extreme_->approximation();

    // If we came from (say) the left, and already found a minimum there, then mark left as explored.
    if (best_minimum_ != extremes_.end())
      last_extreme_->explored(opposite(hdirection_));

    // Update small_step_ (value() returns an absolute value).
    small_step_ = approximation_ptr_->scale().value();
    Dout(dc::notice, "small_step_ set to " << small_step_);

    // Keep track of the best minimum so far; or abort if this minimum isn't better then one found before.
    if (next_extreme_type_ == ExtremeType::minimum)
    {
      Dout(dc::notice, "new extreme = " << *last_extreme_);
      // We were looking for a minimum.
      ASSERT(last_extreme_->is_minimum());
      if (best_minimum_ == extremes_.end() || best_minimum_->cp_sample().Lw() > last_extreme_->cp_sample().Lw())
      {
        best_minimum_ = last_extreme_;
        Dout(dc::notice, "best_minimum_ set to " << best_minimum_->cp_sample() <<
            " and cubic approximation: " << best_minimum_->approximation());
      }
      if (last_extreme_ != best_minimum_)
      {
        // The new minimum isn't better than what we found already. Stop going into this direction.
        Dout(dc::notice, "The new minimum (at " << last_extreme_->cp_sample() << ") isn't better than what we found already. "
            "Stop going into the direction " << hdirection_ << ".");
        state_ = IterationState::abort_hdirection;
        return false;
      }
    }

    // After finding a maximum we want to find a minimum and visa versa. Change next_extreme_type_.
    next_extreme_type_ = opposite(next_extreme_type_);
    Dout(dc::notice, "next_extreme_type_ is toggled to " << next_extreme_type_ << ".");
    // Invalidate what was returned by find_extreme.
    last_region_ = Region::invalid;
    Dout(dc::notice, "Invalidated last_region_ (set to Region::invalid)");
  }

  // Find all samples that are within the scale range of the current approximation.
  Scale const& scale = approximation_ptr_->scale();
  double const scale_value = scale.value();
  // Should only store positive values. Zero would mean "not initialized", aka - !is_valid().
  ASSERT(!scale.negligible(scale_value));

  double const critical_point_w = scale.critical_point_w();

  std::array<Sample const*, 5> usable_samples;
  int number_of_usable_samples = 0;

  auto add_usable_sample = [=, &usable_samples, &number_of_usable_samples](Sample const* sample){
    double dist = std::abs(sample->w() - critical_point_w);
    // If the current sample is too close or too far away from the critical point, then skip this sample.
    if (dist >= 0.001 * scale_value && dist <= 1.1 * scale_value)
      usable_samples[number_of_usable_samples++] = sample;
  };

  // The approximation should have two samples, and the current() one should be the local extreme.
  // Get the other sample that was used for the approximation.
  Sample const* approximation_sample = &approximation_ptr_->prev();       // This sample might or might not be in the history.
  add_usable_sample(approximation_sample);

  // Run over all samples except the local extreme.
  int local_extreme_index = state_ == IterationState::local_extreme ? 0 : 1;
  for (int i = 1 - local_extreme_index; i < history_.relevant_samples() && number_of_usable_samples < usable_samples.size(); ++i)
  {
    if (i == local_extreme_index)
      continue;
    Sample const* sample = &history_.prev(i);
    if (sample == approximation_sample)
      continue;
    add_usable_sample(sample);
  }

  Dout(dc::notice, "Number of samples within scale range: " << number_of_usable_samples);

  // If we do not already have at least three relevant samples then get another one.
  if (number_of_usable_samples < 2)
  {
    // After already asking for one extra sample we should have enough:
    // the extra sample plus the other (non-extreme) sample from the approximation.
    ASSERT(state_ != IterationState::need_extra_sample);

    // Request an extra sample on the other side of the critical point then that we already have.
    HorizontalDirection direction = last_extreme_->approximation().prev_to_current();
    w += last_extreme_->approximation().scale().step(direction);
    expected_Lw_ = last_extreme_->approximation().at(w);
    Debug(set_algorithm_str(w, "extra sample (keep going (no zeroes))"));

    if (hdirection_ == HorizontalDirection::undecided)
    {
      hdirection_ = direction;
      Dout(dc::notice, "Initialized hdirection_ to " << hdirection_ << ".");
    }

    // Can this ever fail?
    ASSERT(hdirection_ == direction);

    // The extra sample should be added on the other side of the critical point.
    ASSERT(number_of_usable_samples == 0 || (w > critical_point_w) != (usable_samples[0]->w() > critical_point_w));

    state_ = IterationState::need_extra_sample;
    return true;
  }

  Sample const& w2 = history_.prev(local_extreme_index);
  double const w2_1 = w2.w();

  // Brute force find the two samples that, together with the local extreme sample, have the largest spread.
  int i0 = 0;
  int i1 = 1;
  double best_spread = 0.0;
  if (number_of_usable_samples > 2)
  {
    for (int t0 = 0; t0 < number_of_usable_samples - 1; ++t0)
      for (int t1 = t0 + 1; t1 < number_of_usable_samples; ++t1)
      {
        double w0 = usable_samples[t0]->w();
        double w1 = usable_samples[t1]->w();
        double spread = utils::square(w0 - w1) + utils::square(w0 - w2_1) + utils::square(w1 - w2_1);
        if (spread > best_spread)
        {
          best_spread = spread;
          i0 = t0;
          i1 = t1;
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
  auto quotient = derivative.long_division(w2_1, remainder);
  Dout(dc::notice, "quotient = " << quotient << " with remainder " << remainder);

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{quotient_event, quotient});
#endif

  std::array<double, 2> zeroes;
  int number_of_zeroes = quotient.get_zeroes(zeroes);

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
    Dout(dc::notice, "with zeroes " << zeroes[0] << " and " << zeroes[1]);
  else if (number_of_zeroes == 1)
    Dout(dc::notice, "with one zero at " << zeroes[0]);
  else
    Dout(dc::notice, "with no zeroes!");

  if (number_of_zeroes > 0)
  {
    std::array<double, 2> expected_Lw;
    int best_zero = 0;
    if (number_of_zeroes == 2)
    {
      if (w2_1 > zeroes[0] == w2_1 > zeroes[1])
      {
        // If both zeroes are on the same side of the found local extreme, pick the nearest one.
        best_zero = std::abs(w2_1 - zeroes[0]) < std::abs(w2_1 - zeroes[1]) ? 0 : 1;
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
      // There is only one usable zero.
      expected_Lw[0] = fourth_degree_approximation(zeroes[0]);
    }
    w = zeroes[best_zero];
    expected_Lw_ = expected_Lw[best_zero];
    Debug(set_algorithm_str(w, "best zero"));
    reset_history();
    Dout(dc::notice(hdirection_ == HorizontalDirection::undecided && number_of_zeroes == 2),
        "Best zero was " << zeroes[best_zero] << " with A(" << zeroes[best_zero] << ") = " <<
        fourth_degree_approximation(zeroes[best_zero]) <<
        " (the other has value " << fourth_degree_approximation(zeroes[1 - best_zero]) << ")");
  }
  else if (state_ == IterationState::need_extra_sample)
  {
    // This should never assert: we already set hdirection_ when requesting the extra sample.
    ASSERT(hdirection_ != HorizontalDirection::undecided);

    // The next call to find_extreme will use this local extreme and the w value that we return
    // as the two samples for the cubic. Then we are only interested in extremes in the same
    // direction as hdirection_.
    hrestriction_ = static_cast<Restriction>(hdirection_);
    Dout(dc::notice, "hrestriction_ is now " << hrestriction_);

    // Remember in which direction we travelled from this extreme.
    last_extreme_->explored(hdirection_);

    // Returning false here signals that we did not find any zero.
    return false;
  }
  else if (hdirection_ == HorizontalDirection::undecided)
  {
    // Keep going in the same direction.
    w += last_extreme_->approximation().scale().step(hrestriction_);
    expected_Lw_ = last_extreme_->approximation().at(w);
    Debug(set_algorithm_str(w, "past extreme (no zeroes)"));
  }
  else
  {
    // Keep going in the same hdirection.
    w += last_extreme_->approximation().scale().step(hdirection_);
    expected_Lw_ = last_extreme_->approximation().at(w);
    Debug(set_algorithm_str(w, "keep going (no zeroes)"));
  }

  if (hdirection_ == HorizontalDirection::undecided)
  {
    // Now that the decision on which hdirection_ we explore is taken, store that decision.
    hdirection_ = w < w2_1 ? HorizontalDirection::left : HorizontalDirection::right;
    Dout(dc::notice, "Initialized hdirection_ to " << hdirection_ << ".");
  }

  // The next call to find_extreme will use this local extreme and the w value that we return
  // as the two samples for the cubic. Then we are only interested in extremes in the same
  // direction as hdirection_.
  hrestriction_ = static_cast<Restriction>(hdirection_);
  Dout(dc::notice, "hrestriction_ is now " << hrestriction_);

  // Remember in which direction we travelled from this extreme.
  last_extreme_->explored(hdirection_);

  // The local extreme was handled.
  state_ = IterationState::done;
  return true;
}

double Algorithm::update_approximation(bool current_is_replacement)
{
  DoutEntering(dc::notice, "Algorithm::update_approximation(" << std::boolalpha << current_is_replacement << ")");

  using namespace gradient_descent;

  bool const first_call = last_region_ == Region::unknown;
  double new_w CWDEBUG_ONLY(= uninitialized_magic);

#ifdef CWDEBUG
  math::CubicPolynomial old_cubic = approximation_ptr_->cubic();
#endif
  ScaleUpdate result;
  if (approximation_ptr_ == &current_approximation_)
  {
    approximation_ptr_->add(&history_.current(), current_is_replacement, next_extreme_type_, state_ == IterationState::back_tracking);

    if (first_call && approximation_ptr_->number_of_relevant_samples() == 2)
    {
      ASSERT(next_extreme_type_ == ExtremeType::unknown);
      new_w = approximation_ptr_->find_extreme(last_region_, next_extreme_type_, hrestriction_);
      // This is expected to be set the first time.
      ASSERT(last_region_ != Region::unknown);
    }

    result = approximation_ptr_->update_scale(current_is_replacement, next_extreme_type_);
  }
  else
  {
    // Does this ever happen?
    ASSERT(state_ != IterationState::back_tracking);
    result = approximation_ptr_->update_local_extreme_scale(history_.current());

    // If this fails then apparently we even have to call find_extreme in this case!
    ASSERT(!(first_call && approximation_ptr_->number_of_relevant_samples() != 1));
  }

  if (result == ScaleUpdate::disconnected)
  {
    // disconnected can only be returned by update_local_extreme_scale.
    ASSERT(approximation_ptr_ != &current_approximation_);

    Sample const* nearest_sample = approximation_ptr_->scale().get_nearest_sample(&history_.current());
    reset_history();    // This sets approximation_ptr_ = &current_approximation_, therefore call add().

    // Start the new approximation with nearest_sample + current.
    approximation_ptr_->add(nearest_sample, false, next_extreme_type_, false);
    approximation_ptr_->add(&history_.current(), false, next_extreme_type_, false);
    result = approximation_ptr_->update_scale(false, next_extreme_type_);
    ASSERT(result == ScaleUpdate::initialized);
  }
#ifdef CWDEBUG
  else
    event_server_.trigger(AlgorithmEventType{scale_draw_event, result, approximation_ptr_->scale(), old_cubic});

  Dout(dc::notice, "approximation = " << *approximation_ptr_ <<
      " (" << utils::print_using(*approximation_ptr_, &Approximation::print_based_on) << ")");

  event_server_.trigger(AlgorithmEventType{cubic_polynomial_event, approximation_ptr_->cubic()});
#endif

  return new_w;
}

void Algorithm::handle_single_sample(Weight& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_single_sample(" << w << ")");

  double step;
#ifdef CWDEBUG
  char const* algorithm_str;
#endif

  if (small_step_ == 0.0)       // Not defined yet?
  {
    step = learning_rate_ * -history_.current().dLdw();

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
    // small_step_ should only be set once hdirection_ has been set.
    ASSERT(hdirection_ != HorizontalDirection::undecided);
    step = static_cast<int>(hdirection_) * small_step_;
#ifdef CWDEBUG
    algorithm_str = "small step";
#endif
  }

  // This step could still be too small.
  if (approximation_ptr_->scale().negligible(step))
  {
    // This wouldn't be good because then the new sample will replace
    // the current one and we'd still have just one sample.
    step = approximation_ptr_->scale().make_significant(step);
#ifdef CWDEBUG
    algorithm_str = "avoiding replacement";
#endif
  }
  w += step;
  expected_Lw_ = approximation_ptr_->at(w);
  Debug(set_algorithm_str(w, algorithm_str));
  state_ = IterationState::done;
}

void Algorithm::handle_approximation(Weight& w, bool first_call, double new_w)
{
  DoutEntering(dc::notice, "Algorithm::handle_approximation(" << w << ", " << std::boolalpha << first_call <<
      ", " << new_w << ")");

  Approximation& approximation(*approximation_ptr_);

  // Is this the first call, or after a reset?
  if (first_call)       // In this case find_extreme is already called and new_w has already been initialized.
  {
    // Call find_extreme before getting here.
    ASSERT(last_region_ != Region::invalid);

    // Sort the two samples that we have, so that the one in the direction of hdirection_ appears "last".
    if (last_region_ != Region::inbetween)
      approximation_ptr_->set_current_index(last_region_);

    Debug(expected_Lw_ = 1234567.89);         // Should not be used.
    if (next_extreme_type_ == ExtremeType::unknown)
    {
      // The third degree polynomial fit does not have local extremes.
      // In this case last_region_ is set to point in the direction where we descent.
      ASSERT(last_region_ != Region::inbetween);
      w += static_cast<int>(last_region_) * approximation_ptr_->scale().step(last_region_);
      Debug(set_algorithm_str(w, "downhill"));
      next_extreme_type_ = ExtremeType::minimum;
      Dout(dc::notice, "Setting next_extreme_type_ to " << next_extreme_type_ << " because we're going downhill.");
      state_ = IterationState::done;
      return;
    }
    else
    {
      // See comment above; find_extreme should have set new_w in this case.
      ASSERT(new_w != uninitialized_magic);
      Debug(set_algorithm_str(new_w, "initial find_extreme jump"));
      state_ = IterationState::extreme_jump;
    }

    // The extreme type set by find_extreme (or the minimum because we're going "downhill") must correspond with the scale type!
    ASSERT((next_extreme_type_ == ExtremeType::minimum) == (approximation.scale().type() == CriticalPointType::minimum) &&
           (next_extreme_type_ == ExtremeType::maximum) == (approximation.scale().type() == CriticalPointType::maximum));
  }
  else
  {
#if CW_DEBUG
    Region prev_region = last_region_;
#endif
    ExtremeType extreme_type = next_extreme_type_;
    new_w = approximation_ptr_->find_extreme(last_region_, extreme_type, hrestriction_);
    // We really shouldn't change our mind about the direction.
    ASSERT(prev_region == Region::invalid || last_region_ == Region::inbetween || last_region_ == prev_region);

    if (extreme_type != ExtremeType::unknown)
    {
      if (next_extreme_type_ == ExtremeType::unknown)
        next_extreme_type_ = extreme_type;
      expected_Lw_ = approximation_ptr_->at(new_w);
      Debug(set_algorithm_str(new_w, "find_extreme jump"));
      // If the extreme is in between the two samples used for the approximation then we jumped too far.
      state_ = last_region_ == Region::inbetween ? IterationState::back_tracking : IterationState::extreme_jump;
    }
    else
    {
      // This means there isn't an extreme of the type that we want. Keep going in the same direction.
      // We want to "keep going" in the same direction. But what to do if hdirection_ is undecided?
      ASSERT(hdirection_ != HorizontalDirection::undecided);
      w += approximation_ptr_->scale().step(hdirection_);
      expected_Lw_ = approximation_ptr_->at(w);
      Debug(set_algorithm_str(w, "keep going"));
      state_ = IterationState::check_energy;
      return;
    }
  }

  // The above code must always set a new w value (or use the one passed).
  ASSERT(new_w != uninitialized_magic);
  double step = w - new_w;
  w = new_w;
  expected_Lw_ = approximation_ptr_->at(w);

  double abs_step = std::abs(step);
  Dout(dc::notice, "abs_step = " << abs_step << " (between " <<
      (history_.total_number_of_samples() - 1) << " and " << history_.total_number_of_samples() << " (to be added))");

  // Did we reach the (local) extreme?
  if (abs_step < (next_extreme_type_ == ExtremeType::minimum ? 0.01 : 0.05) * approximation.scale().value())
  {
    Dout(dc::notice, (next_extreme_type_ == ExtremeType::minimum ? "Minimum" : "Maximum") << " reached: " << abs_step <<
        " < " << (next_extreme_type_ == ExtremeType::minimum ? 0.01 : 0.05) << " * " << approximation.scale().value() <<
        "[" << approximation.scale() << "]");
    state_ = IterationState::local_extreme;
  }
}

bool Algorithm::handle_abort_hdirection(Weight& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_abort_hdirection(" << w << ")");

  using namespace gradient_descent;

  // Has a minimum been found at all?
  if (best_minimum_ == extremes_.end())
  {
    Dout(dc::notice, "Aborting going " << hdirection_ << " and terminating search because no minimum has has been found at all!");
    return false;
  }

  // Jump back to the best minimum and continue in the opposite hdirection.
  Dout(dc::notice, "Aborting exploring " << hdirection_ << " of the minimum at " << best_minimum_->cp_sample() << ".");

  // Was this minimum already explored in both directions?
  if (best_minimum_->done())
    return false;

  ASSERT(hdirection_ != HorizontalDirection::undecided);

  // Change hdirection.
  hdirection_ = opposite(hdirection_);
  Dout(dc::notice, "Changed horizontal direction to " << hdirection_);

  // Now going in a different direction, from a local minimum (which can happen if we
  // didn't explore the next minimum into that direction yet) we must update hrestriction_.
  hrestriction_ = static_cast<Restriction>(hdirection_);
  Dout(dc::notice, "hrestriction_ is now " << hrestriction_);

  // See if we already know a neighbor (maximum) into this direction.
  // This can be the case if we didn't explore into hdirection yet from that maximum.
  LocalExtreme* neighbor = best_minimum_->neighbor(hrestriction_);
  LocalExtreme* new_local_extreme = neighbor ? neighbor : &*best_minimum_;

  // Restore the current sample and scale to the values belonging to this new local extreme.
  w = new_local_extreme->cp_sample().w();
  Dout(dc::notice, "Restored w to best minimum at " << w);

  // Restore small_step_ to what is was.
  small_step_ = new_local_extreme->approximation().scale().value();

  // Do a step with a size equal to the scale in this minimum: we expect it not to drastically change before that point.
  w += static_cast<int>(hdirection_) * small_step_;
  expected_Lw_ = new_local_extreme->approximation().at(w);

  if (neighbor)
  {
    Debug(set_algorithm_str(w, "neighbor of best minimum, opposite direction"));
    next_extreme_type_ = ExtremeType::minimum;
    neighbor->explored(hdirection_);
  }
  else
  {
    Debug(set_algorithm_str(w, "best minimum, opposite direction"));
    next_extreme_type_ = ExtremeType::maximum;
  }

  best_minimum_->explored(hdirection_);
  Dout(dc::notice, "next_extreme_type_ is set to " << next_extreme_type_ << " because we just jumped to a minimum.");

  // Restore the energy to what it was when this minimum was stored.
  energy_.set(new_local_extreme->energy(), new_local_extreme->cp_sample().Lw());

  reset_history();
  last_region_ = Region::invalid;
  Dout(dc::notice, "Invalidated last_region_ (set to Region::invalid)");
  approximation_ptr_->add(&new_local_extreme->cp_sample(), false, next_extreme_type_, false);

  // w was successfully updated.
  Dout(dc::notice, "Returning: " << std::setprecision(std::numeric_limits<double>::max_digits10) << w);
  state_ = IterationState::done;
  return true;
}

#ifdef CWDEBUG
void DifferenceEventData::print_on(std::ostream& os) const
{
  os << "DifferenceEventData:{w:" << w_ << ", expected_Lw:" << expected_Lw_ << ", Lw:" << Lw_ << "}";
}

void PolynomialEventData::print_on(std::ostream& os) const
{
  os << "PolynomialEventData:{" << polynomial_ << "}";
}

void QuadraticPolynomialEventData::print_on(std::ostream& os) const
{
  os << "QuadraticPolynomialEventData:{" << quadratic_polynomial_ << "}";
}

void CubicPolynomialEventData::print_on(std::ostream& os) const
{
  os << "CubicPolynomialEventData:{" << cubic_polynomial_ << "}";
}

void AlgorithmEventData::print_on(std::ostream& os) const
{
  if (std::holds_alternative<ResetEventData>(event_data_))
    std::get<ResetEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<DifferenceEventData>(event_data_))
    std::get<DifferenceEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<FourthDegreeApproximationEventData>(event_data_))
    std::get<FourthDegreeApproximationEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<DerivativeEventData>(event_data_))
    std::get<DerivativeEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<QuotientEventData>(event_data_))
    std::get<QuotientEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<QuadraticPolynomialEventData>(event_data_))
    std::get<QuadraticPolynomialEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<KineticEnergyEventData>(event_data_))
    std::get<KineticEnergyEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<ScaleDrawEventData>(event_data_))
    std::get<ScaleDrawEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<ScaleEraseEventData>(event_data_))
    std::get<ScaleEraseEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<HistoryAddEventData>(event_data_))
    std::get<HistoryAddEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<CubicPolynomialEventData>(event_data_))
    std::get<CubicPolynomialEventData>(event_data_).print_on(os);
  else
    // Missing implementation.
    ASSERT(false);
}

//static
std::string Algorithm::to_string(IterationState state)
{
  switch (state)
  {
    AI_CASE_RETURN(IterationState::done);
    AI_CASE_RETURN(IterationState::check_energy);
    AI_CASE_RETURN(IterationState::extreme_jump);
    AI_CASE_RETURN(IterationState::back_tracking);
    AI_CASE_RETURN(IterationState::local_extreme);
    AI_CASE_RETURN(IterationState::need_extra_sample);
    AI_CASE_RETURN(IterationState::abort_hdirection);
  }
  AI_NEVER_REACHED
}

std::ostream& operator<<(std::ostream& os, Algorithm::IterationState state)
{
  os << Algorithm::to_string(state);
  return os;
}

#endif

} // namespace gradient_descent
