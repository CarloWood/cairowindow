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
  Dout(dc::notice, "hdirection_ = " << hdirection_ << "; vdirection_ = " << vdirection_);

  using namespace gradient_descent;

#if 0 // OLD CODE
  double gamma_based_w;
  if (state_ == IterationState::vertex_jump)
  {
    Approximation& approximation(*approximation_ptr_);
    // How else could we have made a parabolic approximation?
    ASSERT(approximation.number_of_relevant_samples() > 0);
    // Then state should have been IterationState::local_extreme.
    ASSERT(!approximation.is_extreme());

    // See https://math.stackexchange.com/questions/4923841/ plus answer.

    double h = Lw - expected_Lw_;
    double c = approximation.parabola()[2];
    double w0 = approximation.current().w();
    double w1 = w;

    double gamma = 6.0 * h / (c * utils::square(w0 - w1));

    Dout(dc::notice, "h = " << h << "; gamma = " << gamma);

    // Did we substantially miss our target?
    if (gamma > -1.0 && std::abs(gamma) > 0.1)
    {
      // Approximate a better w based on the minimum of a third degree polynomial fit, using the new sample.
      gamma_based_w = w0 + 2.0 * (std::sqrt(1.0 + gamma) - 1.0) / gamma * (w1 - w0);
      //expected_Lw_ = ; Not set... is currently ignored anyway.
      if (gamma > 1.207)
      {
        // Too far off, don't even add w to the history.
        Dout(dc::notice, w << " --> " << history_.total_number_of_samples() << ": " <<
            gamma_based_w << " [gamma based] [expected_Lw: " << expected_Lw_ << "]");
        w = gamma_based_w;
        state_ = IterationState::done;
        Dout(dc::notice, "Returning: " << std::setprecision(std::numeric_limits<double>::max_digits10) << w);
        return true;
      }
      Dout(dc::notice, w << " --> " << (history_.total_number_of_samples() + 1) << ": " <<
          gamma_based_w << " [gamma based] [expected_Lw: " << expected_Lw_ << "]");
      state_ = IterationState::gamma_based;
    }
  }
#endif

  // If the new sample (w) is too close to the previous sample (Scale::negligible returns true)
  // then the new sample replaces the previous sample, current_is_replacement is set to true
  // and no new sample is added to the history.
  //
  // Otherwise, the new sample is added to the history and current_is_replacement is set to false.
  bool current_is_replacement;
  history_.add(w, Lw, dLdw, approximation_ptr_->parabola_scale(), current_is_replacement);

  // This function should never return a value whose difference with the previous sample is negligible
  // if there is only a single relevant sample in the history.
  ASSERT(!current_is_replacement || history_.relevant_samples() > 1);

#ifdef CWDEBUG
  // Erase all previous curves (if they exist).
  event_server_.trigger(AlgorithmEventType{reset_event});
  event_server_.trigger(AlgorithmEventType{difference_event, w, expected_Lw_, Lw});
#endif

  // Update kinetic energy. Returns false if too much energy was used.
  if (!update_energy())
    return handle_abort_hdirection(w);

  // Handle the case where this sample is a local extreme.
  if (state_ == IterationState::local_extreme)
  {
    // This state is set when locally the curve looks like a parabola that we're trying to find the
    // extreme of, and the last sample is that extreme (that is, the derivative is now close to zero).
    // Returns false if this exreme is a minimum but isn't better than the previously found best minimum.
    if (!handle_local_extreme(w))
      return handle_abort_hdirection(w);
    Dout(dc::notice, "Returning: " << std::setprecision(std::numeric_limits<double>::max_digits10) << w);
    return true;        // w was successfully updated by handle_local_extreme.
  }

  // Create/update a parabolic approximation from this and the previous sample (or a line if this is the first sample).
  update_approximation(current_is_replacement);
  Approximation& approximation(*approximation_ptr_);

  if (approximation.number_of_relevant_samples() == 1)
  {
    handle_single_sample(w);
    Dout(dc::notice, "Returning: " << std::setprecision(std::numeric_limits<double>::max_digits10) << w);
    return true;
  }

#if 0 // OLD CODE
  if (state_ == IterationState::gamma_based)
  {
    w = gamma_based_w;
    state_ = IterationState::done;
  }
  else
#endif

  if (!handle_approximation(w))
    return false;
#if 0 // OLD CODE
  if (state_ == IterationState::clamped)
  {
    Dout(dc::notice, "Clamped by sample " << history_[clamped_history_index_]);
    // The vertex of the current approximation doesn't make sense, aka the approximation doesn't make sense.
    // We need to find a reasonable jump point based on the current sample and the one given by clamped_history_index_.
    //FIXME: make this 0.5
    w = 0.52 * (w + history_[clamped_history_index_].w());
    state_ = IterationState::done;
  }
  else
#endif
  if (state_ == IterationState::local_extreme)
  {
    // Do not return a value that would replace the current sample.
    if (approximation_ptr_->parabola_scale().negligible(w - history_.current().w()))
    {
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
    // This state is set when locally the curve looks like a parabola with its minimum
    // on the side that we just came from: in this case we move away from the vertex,
    // losing energy. Since we're going uphill, we need to check if we didn't go too
    // high for the kinetic energy that we have, in which case the exploration of this
    // direction is aborted.

    // If the horizontal direction is still unknown, then we should always go towards the extreme
    // of the local parabolic approximation: in fact we would have jumped there and would not
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
  // Reset the parabolic approximation.
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

  // Store it as an extreme.
  extremes_type::iterator new_extreme =
    extremes_.emplace(hdirection_ == HorizontalDirection::right ? extremes_.end() : extremes_.begin(),
        history_.current(), *approximation_ptr_, energy_.energy());

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{scale_erase_event});
#endif
  // Switch approximation_ptr to the parabolic approximation stored in this extreme:
  // we need to keep updating it when new samples are added that match the same parabolic.
  approximation_ptr_ = &new_extreme->approximation();

  // If we came from (say) the left, and already found a minimum there, then mark left as explored.
  if (best_minimum_ != extremes_.end())
    new_extreme->explored(opposite(hdirection_));

  // Update small_step_ (parabola_scale() returns an absolute value).
  small_step_ = approximation_ptr_->parabola_scale();
  Dout(dc::notice, "small_step_ set to " << small_step_);

  // Keep track of the best minimum so far; or abort if this minimum isn't better then one found before.
  if (vdirection_ == VerticalDirection::down)
  {
    // With vdirection_ down, we were looking for a minimum.
    ASSERT(new_extreme->is_minimum());
    if (best_minimum_ == extremes_.end() || best_minimum_->vertex_sample().Lw() > new_extreme->vertex_sample().Lw())
    {
      best_minimum_ = new_extreme;
      Dout(dc::notice, "best_minimum_ set to " << best_minimum_->vertex_sample() <<
          " and parabolic approximation: " << best_minimum_->approximation());
    }
    if (new_extreme != best_minimum_)
    {
      // The new minimum isn't better than what we found already. Stop going into this direction.
      Dout(dc::notice, "The new minimum (at " << new_extreme->vertex_sample() << ") isn't better than what we found already. "
          "Stop going into the direction " << hdirection_ << ".");
      state_ = IterationState::abort_hdirection;
      return false;
    }
  }

  // After finding a maximum we want to find a minimum and visa versa. Change vdirection_.
  vdirection_ = opposite(vdirection_);
  Dout(dc::notice, "vdirection_ is toggled to " << vdirection_ << ".");

  Sample const& w2 = history_.current();
  double w2_1 = w2.w();

  // Find all samples that are within the scale range of the current approximation.
  double const scale = approximation_ptr_->parabola_scale();
  std::array<int, 5> usable_samples;
  int number_of_usable_samples = 0;
  for (int i = 1; i < history_.relevant_samples() && number_of_usable_samples < usable_samples.size(); ++i)
  {
    Sample const& sample = history_.prev(i);
    double dist = std::abs(sample.w() - w2_1);
    // If the current sample is too close or too far away from the vertex, then skip this sample.
    if (dist < 0.001 * scale || dist > 1.1 * scale)
      continue;

    usable_samples[number_of_usable_samples++] = i;
  }
  Dout(dc::notice, "Number of samples within scale range: " << number_of_usable_samples);
  // If not enough samples we need to get another one! (to be implemented)
  ASSERT(number_of_usable_samples >= 2);
  // Brute force find the two samples that, together with the current sample, have the largest spread.
  int i0 = 0;
  int i1 = 1;
  double best_spread = 0.0;
  if (number_of_usable_samples > 2)
  {
    for (int t0 = 0; t0 < number_of_usable_samples - 1; ++t0)
      for (int t1 = t0 + 1; t1 < number_of_usable_samples; ++t1)
      {
        double w0 = history_.prev(usable_samples[t0]).w();
        double w1 = history_.prev(usable_samples[t1]).w();
        double spread = utils::square(w0 - w1) + utils::square(w0 - w2_1) + utils::square(w1 - w2_1);
        if (spread > best_spread)
        {
          best_spread = spread;
          i0 = t0;
          i1 = t1;
        }
      }
  }

  Sample const& w1 = history_.prev(usable_samples[i0]);
  Sample const& w0 = history_.prev(usable_samples[i1]);

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

  // If we have (at least) three points, the approximation is a fourth degree parabola:
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

  // It is possible that the zeroes are no usable because they are on the wrong side.
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
    for (int zero = 0; zero < number_of_zeroes; ++zero)
      expected_Lw[zero] = fourth_degree_approximation(zeroes[zero]);
    int best_zero = (number_of_zeroes == 2 &&
        (hdirection_ == HorizontalDirection::right ||
         (hdirection_ == HorizontalDirection::undecided &&
          expected_Lw[1] < expected_Lw[0]))) ? 1 : 0;
    w = zeroes[best_zero];
    expected_Lw_ = expected_Lw[best_zero];
    Debug(set_algorithm_str(w, "best zero"));
    reset_history();
    Dout(dc::notice(hdirection_ == HorizontalDirection::undecided && number_of_zeroes == 2),
        "Best zero was " << zeroes[best_zero] << " with A(" << zeroes[best_zero] << ") = " <<
        fourth_degree_approximation(zeroes[best_zero]) <<
        " (the other has value " << fourth_degree_approximation(zeroes[1 - best_zero]) << ")");
  }
  else if (hdirection_ == HorizontalDirection::undecided)
  {
    // The "scale" of the (current) parabolic approximation is set to 'edge sample' minus x-coordinate of the vertex (v_x),
    // where 'edge sample' is the sample that is "part of" the parabolic approximation that is the furthest away
    // from the vertex. "Part of" here means that it deviates vertically less than 10% of the vertical distance to
    // the vertex. In that sense we can consider that we "came from" the direction of that edge sample.
    // For example, if the edge sample had x-coordinate w = w_E and we came from the right (w_E > v_x) then
    // scale = w_E - v_x > 0. Subtracting a positive scale thus means we continue in the direction towards to the left.
    // If w_E is on the left of the vertex then scale is negative and subtracting is causes us to continue to the right.
    //
    // Keep going in the same direction.
    w -= new_extreme->approximation().parabola_scale();
    expected_Lw_ = new_extreme->approximation().parabola()(w);
    Debug(set_algorithm_str(w, "past extreme (no zeroes)"));

    // Note that w was already set to the v_x before, but the test below still works
    // because |v_x - history_.current().w()| was determined to be less than 1% of approximation.parabola_scale().
  }
  else
  {
    // Keep going in the same hdirection.
    w += static_cast<int>(hdirection_) * std::abs(new_extreme->approximation().parabola_scale());
    expected_Lw_ = new_extreme->approximation().parabola()(w);
    Debug(set_algorithm_str(w, "keep going (no zeroes)"));
  }

  if (hdirection_ == HorizontalDirection::undecided)
  {

    // Now that the decision on which hdirection_ we explore is taken, store that decision.
    hdirection_ = w - history_.current().w() < 0.0 ? HorizontalDirection::left : HorizontalDirection::right;
    Dout(dc::notice, "Initialized hdirection_ to " << hdirection_ << ".");
  }

  // Remember in which direction we travelled from this extreme.
  new_extreme->explored(hdirection_);

  // The local extreme was handled.
  state_ = IterationState::done;
  return true;
}

void Algorithm::update_approximation(bool current_is_replacement)
{
  DoutEntering(dc::notice, "Algorithm::update_approximation(" << std::boolalpha << current_is_replacement << ")");

  using namespace gradient_descent;

  math::QuadraticPolynomial old_parabola = approximation_ptr_->parabola();
  ScaleUpdate result;
  if (approximation_ptr_ == &current_approximation_)
    result = approximation_ptr_->add(&history_.current(), current_is_replacement);
  else
    result = approximation_ptr_->update_scale(history_.current());

  if (result == ScaleUpdate::disconnected)
  {
#if 0
    reset_history();                                                    // This sets approximation_ptr_ = &current_approximation_.
    result = approximation_ptr_->add(&history_.current(), false);       // Therefore call add().
#endif
  }
#ifdef CWDEBUG
  else
  {
    event_server_.trigger(AlgorithmEventType{scale_draw_event, result,
        approximation_ptr_->parabola_scale().parabola().vertex_x(),
        approximation_ptr_->parabola_scale().edge_sample_w(),
        old_parabola});
  }
#endif

  Dout(dc::notice, "approximation = " << *approximation_ptr_ <<
      " (" << utils::print_using(*approximation_ptr_, &Approximation::print_based_on) << ")");

#ifdef CWDEBUG
  event_server_.trigger(AlgorithmEventType{quadratic_polynomial_event, approximation_ptr_->parabola()});
  if (approximation_ptr_->number_of_relevant_samples() == 2)
    event_server_.trigger(AlgorithmEventType{cubic_polynomial_event, approximation_ptr_->cubic()});
#endif
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
  if (approximation_ptr_->parabola_scale().negligible(step))
  {
    // This wouldn't be good because then the new sample will replace
    // the current one and we'd still have just one sample.
    step = approximation_ptr_->parabola_scale().make_significant(step);
#ifdef CWDEBUG
    algorithm_str = "avoiding replacement";
#endif
  }
  w += step;
  expected_Lw_ = approximation_ptr_->parabola()(w);
  Debug(set_algorithm_str(w, algorithm_str));
  state_ = IterationState::done;
}

bool Algorithm::handle_approximation(Weight& w)
{
  DoutEntering(dc::notice, "Algorithm::handle_approximation(" << w << ")");

  Approximation& approximation(*approximation_ptr_);
  double new_w;

  // Is this the first call, or after a reset?
  if (hdirection_ == HorizontalDirection::undecided)
  {
    ASSERT(vdirection_ == VerticalDirection::unknown);
    new_w = approximation_ptr_->find_extreme(hdirection_, vdirection_);

    // find_extreme never returns undecided.
    ASSERT(hdirection_ != HorizontalDirection::undecided);

    // Sort the two samples that we have, so that the one in the direction of hdirection_ appears "last".
    if (hdirection_ != HorizontalDirection::inbetween)
      approximation_ptr_->set_current_index(hdirection_);

    if (vdirection_ == VerticalDirection::unknown)
    {
      // The third degree polynomial fit does not have local extremes.
      // In this case hdirection_ is set to point in the direction where we descent.
      w += static_cast<int>(hdirection_) * std::abs(approximation_ptr_->parabola_scale());
      Debug(set_algorithm_str(w, "downhill"));
      vdirection_ = VerticalDirection::down;
      Dout(dc::notice, "Setting vdirection_ to " << vdirection_ << " because we're going downhill.");
      state_ = IterationState::done;
      return true;
    }
    else
    {
      Debug(set_algorithm_str(new_w, "find_extreme jump"));
    }
  }
  else if (hdirection_ != HorizontalDirection::inbetween)
  {
#if CW_DEBUG
    HorizontalDirection prev_hdirection = hdirection_;
#endif
    new_w = approximation_ptr_->find_extreme(hdirection_, vdirection_);
    // We really shouldn't change our mind about the direction.
    ASSERT(hdirection_ == HorizontalDirection::inbetween || hdirection_ == prev_hdirection);

    // What to do in this case?
    ASSERT(vdirection_ != VerticalDirection::unknown);
  }
  else
  {
    // Implement.
    Dout(dc::warning, "This code is not implemented!");
    return false;
  }

#if 0 // OLD CODE
  // We have two relevant samples, and thus a parabolic approximation.
  //
  // Case:
  //    A        B        C       D
  // \     /  \     /    ˏ-ˎ     ˏ-ˎ
  // ↑\   /    \   /↑   /   \   /   \
  //   `-´      `-´    /↑    \ /    ↑\
  //
  //                     Jump to the vertex if:
  //      extreme_type   vdirection
  //   A   down             down
  //   B   down             down
  //   C   up               up
  //   D   up               up
  //
  // If the extreme of the current parabolic approximation matches
  // the value of vdirection_ (that the vertical direction that we
  // want to find the next extreme in) then we will jump to that
  // vertex; hdirection_ is not relevant in that case: that is only
  // used when vdirection doesn't match the extreme.

  Sample const& current = approximation.current();
  Sample const& prev = approximation.prev();
  VerticalDirection extreme_type = approximation.parabola_has_maximum() ? VerticalDirection::up : VerticalDirection::down;

  if (vdirection_ != extreme_type)
  {
    double v_x = approximation.parabola().vertex_x();
    // If hdirection matches the parabola; jump to the vertex before adding the scale to w.
    if (hdirection_ == (v_x < w ? HorizontalDirection::left : HorizontalDirection::right))
      w = v_x;
    // There is no extreme in the direction that we're going.
    // Just keep going in the same direction as before.
    double step = static_cast<int>(hdirection_) * std::abs(approximation.parabola_scale());
    w += step;
    expected_Lw_ = approximation.parabola()(w);
    Debug(set_algorithm_str(w, "continue same direction"));
    Dout(dc::notice, (hdirection_ == HorizontalDirection::left ? "Incremented" : "Decremented") << " w with scale (" << std::abs(step) << ")");
    // Abort if the result requires more energy than we have.
    state_ = IterationState::check_energy;
    return;
  }

  // Set w to the value where the derivative of this parabolic approximation is zero.
  Dout(dc::notice, "Setting new_sample to the extreme of parabolic approximation:");
  bool looking_for_maximum = approximation.parabola_has_maximum();
  auto ignore = [&current, looking_for_maximum](Sample const& sample) {
    double w0 = current.w();
    double Lw0 = current.Lw();
    double dLdw0 = current.dLdw();
    double w1 = sample.w();
    double Lw1 = sample.Lw();
    double dLdw1 = sample.dLdw();

    // See clamp_check.cxx.
    double dw = w0 - w1;
    double dw3 = std::pow(dw, 3.0);
    double d = (-2.0 * (Lw0 - Lw1) + dw * (dLdw0 + dLdw1)) / dw3;
    double c = (dLdw0 - dLdw1) / (2.0 * dw) - 1.5 * (w0 + w1) * d;
    double b = (w0 * dLdw1 - w1 * dLdw0) / dw + 3.0 * w0 * w1 * d;

    double D = utils::square(c) - 3.0 * b * d;
    if (D >= 0.0)
    {
      double zero = (-c + (looking_for_maximum ? -1.0 : 1.0) * std::sqrt(D)) / (3.0 * d);
      if (std::min(w0, w1) < zero && zero < std::max(w0, w1))
        return false;
    }

    // There is no minimum/maximum in between w0 and w1, therefore
    // do not clamp on this sample: ignore it.
    return true;
  };
  double new_w = history_.clamp(w, approximation.parabola().vertex_x(), ignore, clamped_history_index_);
  if (!clamped_history_index_.undefined())
  {
    state_ = IterationState::clamped;
    return;
  }
#endif

  double step = w - new_w;
  w = new_w;
  expected_Lw_ = approximation_ptr_->cubic()(w);
  state_ = IterationState::extreme_jump;

  double abs_step = std::abs(step);
  Dout(dc::notice, "abs_step = " << abs_step << " (between " <<
      (history_.total_number_of_samples() - 1) << " and " << history_.total_number_of_samples() << " (to be added))");

  // Did we reach the (local) extreme?
  if (abs_step < (vdirection_ == VerticalDirection::down ? 0.01 : 0.05) * approximation.parabola_scale())
  {
    Dout(dc::notice, (vdirection_ == VerticalDirection::down ? "Minimum" : "Maximum") << " reached: " << abs_step <<
        " < " << (vdirection_ == VerticalDirection::down ? 0.01 : 0.05) << " * " << approximation.parabola_scale());
    // If we do not already have at least three relevant samples, then delay reporting a local extreme
    // until after adding one more sample, making sure it won't replace a previous one.
    if (history_.relevant_samples() < 3)
    {
      // Instead of returning this extreme, do a little overshoot.
      double step = static_cast<int>(hdirection_) * Scale::epsilon;
      w += approximation.parabola_scale().make_significant(step);
      expected_Lw_ = approximation_ptr_->cubic()(w);
      Debug(set_algorithm_str(w, "jump to vertex + small overshoot"));
      return true;
    }
    state_ = IterationState::local_extreme;
  }

  Debug(set_algorithm_str(w, "jump to vertex"));
  return true;
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
  Dout(dc::notice, "Aborting exploring " << hdirection_ << " of the minimum at " << best_minimum_->vertex_sample() << ".");

  // Was this minimum already explored in both directions?
  if (best_minimum_->done())
    return false;

  ASSERT(hdirection_ != HorizontalDirection::undecided);

  // Change hdirection.
  hdirection_ = opposite(hdirection_);
  Dout(dc::notice, "Changed horizontal direction to " << hdirection_);

  // Restore the current sample and scale to the values belonging to this minimum.
  w = best_minimum_->vertex_sample().w();
  Dout(dc::notice, "Restored w to best minimum at " << w);

  // Do a step with a size equal to the scale in this minimum: we expect it not to drastically change before that point.
  w += static_cast<int>(hdirection_) * std::abs(best_minimum_->approximation().parabola_scale());
  expected_Lw_ = best_minimum_->approximation().parabola()(w);
  Debug(set_algorithm_str(w, "best minimum, opposite direction"));
  best_minimum_->explored(hdirection_);
  vdirection_ = VerticalDirection::up;
  Dout(dc::notice, "vdirection_ is set to " << vdirection_ << " because we just jumped to a minimum.");

  // Restore the energy to what it was when this minimum was stored.
  energy_.set(best_minimum_->energy(), best_minimum_->vertex_sample().Lw());

  // Restore small_step_ to what is was.
  small_step_ = best_minimum_->approximation().parabola_scale();

  reset_history();

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
#endif

} // namespace gradient_descent
