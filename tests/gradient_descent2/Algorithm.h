#pragma once

#include "ExtremeChain.h"
#include "IterationState.h"
#include "AlgorithmEventType.h"
#include "KineticEnergy.h"
#include "utils/UniqueID.h"

namespace gradient_descent {

class Algorithm
{
 public:
  // This value is assumed to be much smaller than any change in w over which L(w) changes significantly.
  static constexpr double epsilon = 1e-30;
  static constexpr double negligible_scale_fraction = ExtremeChain::negligible_scale_fraction;
  static constexpr double significant_scale_fraction = ExtremeChain::significant_scale_fraction;

 private:
  double learning_rate_;                // In unit_of(w)^2 / unit_of(L).
  double small_step_{};                 // This will replace learning_rate_ as soon as we have an idea of the scale of changes.
  IterationState state_;
  ExtremeChain chain_;                  // A doubly linked list of SampleNode's, sorted by w value.
  ExtremeType next_extreme_type_;       // The extreme type (minimum or maximum) that we're looking for (next).
  SampleNode::const_iterator left_of_{chain_.end()};    // If not end, then the next extreme (of next_extreme_type_)
                                                        // must found left of this sample.
  SampleNode::const_iterator right_of_{chain_.end()};   // If not end, then the next extreme (of next_extreme_type_)
                                                        // must found right of this sample.
  SampleNode::const_iterator cubic_used_{chain_.end()}; // The node containing the last cubic that was used to jump one of its extremes.
  HorizontalDirection hdirection_;      // The direction relative to FIXME that we want to find the next extreme in.
  KineticEnergy energy_;
  bool have_expected_Lw_{false};        // True if expected_Lw_ was set.
  double expected_Lw_;                  // Whenever w is changed, this is set to what Lw value the approximation is expecting there.

#ifdef CWDEBUG
  double bogus_{-60.0};
  events::Server<AlgorithmEventType> event_server_;
  char const* algorithm_str_;
  utils::UniqueIDContext<int> label_context_;
#endif

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

  bool operator()(double& w, double Lw, double dLdw);
  void handle_single_sample(double& w);
  [[nodiscard]] bool update_energy(double Lw);
  [[nodiscard]] bool handle_abort_hdirection(double& w);
  [[nodiscard]] bool handle_local_extreme(double& w);

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

#ifdef CWDEBUG
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
#endif
};

} // namespace gradient_descent
