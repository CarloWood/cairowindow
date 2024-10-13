#pragma once

#include "ExtremeChain.h"
#include "IterationState.h"
#include "AlgorithmEventType.h"
#include "KineticEnergy.h"
#include "utils/UniqueID.h"
#include "utils/Badge.h"
#ifdef CWDEBUG
#include "debug.h"
#endif

namespace gradient_descent {

class WeightRef;

class Algorithm
{
 public:
  // This value is assumed to be much smaller than any change in w over which L(w) changes significantly.
  static constexpr double epsilon = 1e-30;
  static constexpr double negligible_scale_fraction = ExtremeChain::negligible_scale_fraction;
  static constexpr double significant_scale_fraction = ExtremeChain::significant_scale_fraction;
  static constexpr double humongous_step_scale_factor = 20.0;

 private:
  double learning_rate_;                                // In unit_of(w)^2 / unit_of(L).
  double small_step_{};                                 // This will replace learning_rate_ as soon as we have an idea of the scale of changes.
  IterationState state_;
  ExtremeChain chain_;                                  // A doubly linked list of SampleNode's, sorted by w value.
  ExtremeType next_extreme_type_;                       // The extreme type (minimum or maximum) that we're looking for (next).
  SampleNode::const_iterator left_of_{chain_.end()};    // If not end, then the next extreme (of next_extreme_type_)
                                                        // must found left of this sample.
  SampleNode::const_iterator right_of_{chain_.end()};   // If not end, then the next extreme (of next_extreme_type_)
                                                        // must found right of this sample.
  SampleNode::const_iterator cubic_used_{chain_.end()}; // The node containing the last cubic that was used to jump to one of its extrema.
  HorizontalDirection hdirection_;                      // The direction, relative to the last extreme, that we want to find the next extreme in.
  KineticEnergy energy_;
  bool check_energy_{false};                            // Set to true iff the last probe is a "keep going" step.
  bool have_expected_Lw_{false};                        // True if expected_Lw_ was set.
  double expected_Lw_;                                  // Whenever w is changed, this is set to what Lw value the approximation is expecting there.
  SampleNode::const_iterator last_extreme_cubic_{chain_.end()}; // After calling handle_local_extreme while state_ is not extra_sample,
                                                        // this will point to the cubic that was used to find the local extreme.
                                                        // If that cubic contains two extremes, then it was not marked as local extreme yet.
                                                        // This is also set to point to the local extreme that we jump (back) to, aka before
                                                        // calling handle_local_extreme_jump.
                                                        // Therefore, this can be viewed as being the last local extreme that was visited.
  SampleNode::const_iterator best_minimum_cubic_{chain_.end()}; // Copy of the best last_extreme_cubic_ that was a minimum, so far.
  double best_minimum_energy_;                          // The energy that we had when in the local minimum that is the best minimum.
  bool saw_minimum_ = false;                            // Set to true once we find our first mininum. After it remains true, because
                                                        // we only go to adjacent local extreme, or jump (back) to the best minimum.
                                                        // TODO: reset this to false if we "start over" because we ended up somewhere that
                                                        // is lower than the best minimum so far.

#ifdef CWDEBUG
  events::Server<AlgorithmEventType> event_server_;
  char const* algorithm_str_;
  utils::UniqueIDContext<int> label_context_;
#endif

 private:
  void set_hdirection(HorizontalDirection hdirection)
  {
    hdirection_ = hdirection;
#ifdef CWDEBUG
    Dout(dc::notice, "Set hdirection_ to " << hdirection_ << ".");
    event_server_.trigger(AlgorithmEventType{hdirection_known_event, *last_extreme_cubic_, hdirection_});
#endif
  }

  void initialize_range(double extreme_w);
  void initialize_node(SampleNode::iterator node, SampleNode::const_iterator next COMMA_CWDEBUG_ONLY(bool node_is_last));

  bool wrong_direction(double step) const
  {
    ASSERT(hdirection_ != HorizontalDirection::undecided);
    return (step > 0.0) == (hdirection_ == HorizontalDirection::left);
  }

  void fix_flat(SampleNode::iterator non_const_left_node, SampleNode::iterator non_const_new_node);

 public:
  Algorithm(double learning_rate, double L_max) :
    learning_rate_(learning_rate),
    state_(IterationState::initialization),
    next_extreme_type_(ExtremeType::unknown),
    hdirection_(HorizontalDirection::undecided),
    energy_(L_max COMMA_CWDEBUG_ONLY(event_server_))
  {
    DoutEntering(dc::notice, "Algorithm::Algorithm(" << learning_rate << ", " << L_max << ")");
  }

  bool operator()(WeightRef w, double Lw, double dLdw);
  void handle_single_sample(WeightRef w);
  void move_into_range(WeightRef w);
  [[nodiscard]] bool update_energy(double Lw);
  [[nodiscard]] bool handle_abort_hdirection(WeightRef w);
  [[nodiscard]] bool handle_local_extreme_jump(WeightRef w);
  [[nodiscard]] bool handle_local_extreme(WeightRef w);

  bool success() const
  {
    return state_ == IterationState::success;
  }

  void finish()
  {
    // Get one more sample, then call set_global_minimum.
    state_ = IterationState::finish;
  }

  void set_global_minimum(SampleNode::const_iterator global_minimum)
  {
    cubic_used_ = global_minimum;       // Abuse cubic_used_.
    state_ = IterationState::success;
  }

  Sample const& minimum() const
  {
    ASSERT(state_ == IterationState::success);
    return *cubic_used_;
  }

  // Don't call these from outside of Algorithm. It is only public because WeightRef needs to friend them.
  inline void handle_single_sample_step(WeightRef w, double step, utils::Badge<Algorithm>);
  inline void handle_get_extra_sample(WeightRef w, double requested_w, utils::Badge<Algorithm>);
  inline void handle_divide_last_extreme_cubic(WeightRef w, utils::Badge<Algorithm>);
  inline void extreme_jump_range_edge(WeightRef w, SampleNode::const_iterator range_edge, utils::Badge<Algorithm>);
  inline void handle_extreme_jump(WeightRef w, SampleNode::const_iterator cubic, AnalyzedCubic const& acubic, utils::Badge<Algorithm>);
  inline void handle_last_extreme_plus_step(WeightRef w, utils::Badge<Algorithm>);
  inline void handle_fourth_degree_approximation_jump(WeightRef w, double extreme_w, double extreme_Lw, utils::Badge<Algorithm>);
  inline void handle_finish(WeightRef w, utils::Badge<Algorithm>);

#ifdef CWDEBUG
 public:
  // Accessors for the event servers.
  events::Server<AlgorithmEventType>& event_server()
  {
    return event_server_;
  }

  // Accessors for testsuite.
  double debug_small_step() const { return small_step_; }
  HorizontalDirection debug_hdirection() const { return hdirection_; }
  ExtremeType debug_next_extreme_type() const { return next_extreme_type_; }
  std::string algorithm_str() const { return algorithm_str_; }
  ExtremeChain const& debug_chain() const { return chain_; }
  ExtremeChain& debug_chain() { return chain_; }

  // Manipulators for the testsuite.
  void set_algorithm_str(double new_w, char const* algorithm_str);
  void debug_set_hdirection_next_extreme_type_small_step(HorizontalDirection hdirection, ExtremeType next_extreme_type, double small_step)
  {
    hdirection_ = hdirection;
    next_extreme_type_ = next_extreme_type;
    small_step_ = small_step;
  }

  SampleNode::const_iterator debug_left_of() const { return left_of_; }
  SampleNode::const_iterator debug_right_of() const { return right_of_; }
  SampleNode::const_iterator debug_cubic_used() const { return cubic_used_; }

  bool debug_within_range(double w) const
  {
    return (left_of_ == chain_.end() || w <= left_of_->w()) && (right_of_ == chain_.end() || right_of_->w() <= w);
  }

  void print_range_on(std::ostream& os) const;
#endif
};

class WeightRef
{
 private:
  double& ref_;

  // These need to be able to directly assign to ref_.
  friend void Algorithm::handle_single_sample_step(WeightRef w, double step, utils::Badge<Algorithm>);
  friend void Algorithm::handle_get_extra_sample(WeightRef w, double requested_w, utils::Badge<Algorithm>);
  friend void Algorithm::handle_divide_last_extreme_cubic(WeightRef w, utils::Badge<Algorithm>);
  friend void Algorithm::extreme_jump_range_edge(WeightRef w, SampleNode::const_iterator range_edge, utils::Badge<Algorithm>);
  friend void Algorithm::handle_extreme_jump(WeightRef w, SampleNode::const_iterator cubic, AnalyzedCubic const& acubic,
      utils::Badge<Algorithm>);
  friend void Algorithm::handle_last_extreme_plus_step(WeightRef w, utils::Badge<Algorithm>);
  friend void Algorithm::handle_fourth_degree_approximation_jump(WeightRef w, double extreme_w, double extreme_Lw, utils::Badge<Algorithm>);
  friend void Algorithm::handle_finish(WeightRef w, utils::Badge<Algorithm>);

 public:
  WeightRef(double& w) : ref_(w) { }

  // Read-only accessor.
  operator double() const { return ref_; }

  void find_extreme_jump(SampleNode::const_iterator cubic, SampleNode const& next_node, ExtremeType& next_extreme_type)
  {
    ref_ = cubic->find_extreme(next_node, next_extreme_type);
  }

  void node_jump_step(SampleNode::const_iterator node, double step)
  {
    ref_ = node->w() + step;
  }

 private:
  void extreme_jump(double extreme_w)
  {
    ref_ = extreme_w;
  }

  void extreme_jump_plus_step(SampleNode::const_iterator last_extreme_cubic, double step)
  {
    ref_ = last_extreme_cubic->extreme_w() + step;
  }

#ifdef CWDEBUG
 public:
  // Allow printing without depending on the above.
  friend std::ostream& operator<<(std::ostream& os, WeightRef w)
  {
    return os << w.ref_;
  }
#endif
};

void Algorithm::handle_single_sample_step(WeightRef w, double step, utils::Badge<Algorithm>)
{
  DoutEntering(dc::notice, "Algorithm::handle_single_sample_step(" << w << ", " << step << ")");
  ASSERT(state_ == IterationState::initialization);

  // Never requires a special action-- just update w.ref_.
  w.ref_ += step;
  have_expected_Lw_ = false;
}

void Algorithm::handle_extreme_jump(WeightRef w, SampleNode::const_iterator cubic, AnalyzedCubic const& acubic, utils::Badge<Algorithm>)
{
  DoutEntering(dc::notice, "Algorithm::handle_extreme_jump(" << w << ", " << *cubic << ", " << acubic << ")");
  // acubic must be initialized from cubic.
  ASSERT(acubic.debug_cubic() == &cubic->cubic());

  cubic_used_ = cubic;
  w.extreme_jump(acubic.get_extreme());

#ifdef CWDEBUG
  // Show the cubic that is being used to jump.
  event_server_.trigger(AlgorithmEventType{cubic_polynomial_event, cubic_used_->cubic(), hdirection_});
#endif
  expected_Lw_ = cubic_used_->cubic()(w);
  have_expected_Lw_ = true;
}

void Algorithm::handle_get_extra_sample(WeightRef w, double requested_w, utils::Badge<Algorithm>)
{
  DoutEntering(dc::notice, "Algorithm::handle_get_extra_sample(" << w << ", " << requested_w << ")");
  // This forced assignment is required to be "in range".
  //FIXME: why? Is it really a problem if it is out of range?
//  ASSERT(debug_within_range(requested_w));

  Dout(dc::notice, "Not enough samples to fit a fourth degree polynomial: asking for another samples at w = " << requested_w);
  // Never requires a special action-- just update w.ref_.
  w.ref_ = requested_w;
  have_expected_Lw_ = false;
  state_ = IterationState::extra_sample;
}

void Algorithm::handle_divide_last_extreme_cubic(WeightRef w, utils::Badge<Algorithm>)
{
  DoutEntering(dc::notice, "Algorithm::handle_divide_last_extreme_cubic(" << w << ")");
  math::CubicPolynomial const& cubic = cubic_used_->cubic();
  w.ref_ = cubic.inflection_point();
  expected_Lw_ = cubic[0] + cubic[2] * (2.0 * utils::square(cubic[2] / 3.0) - cubic[3] * cubic[1]) / (3.0 * utils::square(cubic[3]));
  have_expected_Lw_ = true;
//  state_ = IterationState::divide_last_extreme_cubic;
}

void Algorithm::extreme_jump_range_edge(WeightRef w, SampleNode::const_iterator range_edge, utils::Badge<Algorithm>)
{
  DoutEntering(dc::notice, "Algorithm::extreme_jump_range_edge(" << w << ", [" << range_edge->label() << "])");
  // This function should only be used to jump to the flat end of a range.
  ASSERT(range_edge == left_of_ || range_edge == right_of_);
  w.extreme_jump(range_edge->w());
  expected_Lw_ = range_edge->Lw();
  have_expected_Lw_ = true;
}

void Algorithm::handle_last_extreme_plus_step(WeightRef w, utils::Badge<Algorithm>)
{
  DoutEntering(dc::notice, "Algorithm::handle_last_extreme_plus_step(" << w << ")");
  // Keep going in the same hdirection.
  auto& scale = last_extreme_cubic_->scale();
  w.extreme_jump_plus_step(last_extreme_cubic_, scale.step(hdirection_));
  expected_Lw_ = last_extreme_cubic_->cubic()(w);
  have_expected_Lw_ = true;
}

void Algorithm::handle_fourth_degree_approximation_jump(WeightRef w, double extreme_w, double extreme_Lw, utils::Badge<Algorithm>)
{
  DoutEntering(dc::notice, "Algorithm::handle_fourth_degree_approximation_jump(" << w << ", " << extreme_w << ", " << extreme_Lw << ")");
  w.extreme_jump(extreme_w);
  expected_Lw_ = extreme_Lw;
  have_expected_Lw_ = true;
}

void Algorithm::handle_finish(WeightRef w, utils::Badge<Algorithm>)
{
  DoutEntering(dc::notice, "Algorithm::handle_finish(" << w << ")");
  // There is always a minimum :/.
  ASSERT(best_minimum_cubic_ != chain_.end());
  w.extreme_jump(best_minimum_cubic_->extreme_w());
  finish();
}

} // namespace gradient_descent
