#pragma once

#include "Approximation.h"
#include "HorizontalDirection.h"
#include "VerticalDirection.h"
#include "Weight.h"
#include "Sample.h"
#include "Polynomial.h"
#include "PlotHistory.h"
#include "PlotLocalExtreme.h"
#include "KineticEnergy.h"

namespace gradient_descent {

class Algorithm
{
 private:
  double learning_rate_;        // In unit_of(w)^2 / unit_of(L).
  double small_step_{};         // This will replace learning_rate_ as soon as we have an idea of the scale of changes.

#ifdef CWDEBUG
  events::Server<AlgorithmEventType> event_server_;
#endif

  // Remember the (most recent) history of samples.
  PlotHistory history_;
  HistoryIndex clamped_history_index_;
  Approximation current_approximation_;
  Approximation* approximation_ptr_;
  KineticEnergy energy_;
  double expected_Lw_;                  // Whenever w is changed, this is set to what Lw value the approximation is expecting there.

  // hdirection_ is set when we find a local minimum and decide to explore left or right of that.
  // The result is that we'll rather go away from the vertex of the current matching parabolic approximation
  // then towards it, if that doesn't match the current hdirection_.
  HorizontalDirection hdirection_;
  VerticalDirection vdirection_;

  // Using a std::list: pointers to elements may never be invalidated.
  using extremes_type = std::list<PlotLocalExtreme>;
  extremes_type extremes_;
  extremes_type::iterator best_minimum_;
  extremes_type::iterator last_extreme_;

  enum class IterationState
  {
    done,                     // w was already updated.
    check_energy,             // After adding the new sample, abort if the required energy is too large.
    vertex_jump,              // Only add the new sample to the history if we didn't overshoot another extreme dramatically.
    clamped,                  // Attemped to jump to a vertex, but against hdirection, passed a previous sample. Use that last sample instead.
    gamma_based,              // Internal state used to signify that w was already updated and a call to handle_parabolic_approximation
                              // is not longer desired.
    local_extreme,            // After adding the new sample, handle the fact that we found an extreme.
    abort_hdirection          // Stop going in the current hdirection_.
  };

  IterationState state_{IterationState::done};

 public:
  Algorithm(double learning_rate, double L_max,
      cairowindow::plot::Plot& plot, boost::intrusive_ptr<cairowindow::Layer> const& layer,
      cairowindow::draw::PointStyle const& point_style, cairowindow::draw::TextStyle const& label_style) :
    learning_rate_(learning_rate),
    history_(plot, layer, point_style, label_style),
    approximation_ptr_(&current_approximation_),
    energy_(L_max COMMA_CWDEBUG_ONLY(event_server_)),
    best_minimum_(extremes_.end()),
    last_extreme_(extremes_.end()),
    hdirection_(HorizontalDirection::undecided),
    vdirection_(VerticalDirection::down)
  {
  }

  bool operator()(Weight& w, double Lw, double dLdw);
  bool update_energy();
  void reset_history();
  bool handle_local_extreme(Weight& w);
  void update_approximation(bool current_is_replacement);
  void handle_single_sample(Weight& w);
  void handle_parabolic_approximation(Weight& w);
  bool handle_abort_hdirection(Weight& w);

  bool success() const
  {
    return best_minimum_ != extremes_.end();
  }

  Sample const& minimum() const
  {
    // Only call this function when success() returns true.
    ASSERT(success());
    return best_minimum_->vertex_sample();
  }

#ifdef CWDEBUG
  // Accessors for the event servers.
  events::Server<AlgorithmEventType>& event_server()
  {
    return event_server_;
  }
#endif
};

} // namespace gradient_descent