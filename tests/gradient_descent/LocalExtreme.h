#pragma once

#include "Sample.h"
#include "Approximation.h"
#include "HorizontalDirection.h"
#include "debug.h"

namespace gradient_descent {

class LocalExtreme
{
 protected:
  Approximation approximation_;         // The (parabolic) polynomial approximation around this extreme.
  Sample vertex_sample_;                // Sample taken at the vertex of the approximation_;
  double energy_;                       // The maximum height (Lw) that can be reached from here.
  int explored_{0};                     // Bit 1: exploration to the left of this extreme has started.
                                        // Bit 2: same, on the right.
 public:
  LocalExtreme(Sample const& vertex_sample, Approximation const& approximation, double energy) :
    vertex_sample_(vertex_sample), approximation_(approximation), energy_(energy)
  {
    approximation_.set_is_extreme();
  }

  Approximation& approximation()
  {
    return approximation_;
  }

  bool is_minimum() const
  {
    return approximation_.parabola()[2] > 0.0;
  }

  Sample const& vertex_sample() const { return vertex_sample_; }

  double energy() const { return energy_; }

  void explored(HorizontalDirection hdirection)
  {
    DoutEntering(dc::notice, "LocalExtreme::explored(" << hdirection << ") for extreme at " << vertex_sample_ << '.');
    ASSERT(hdirection != HorizontalDirection::undecided);
    int explore_flag = hdirection == HorizontalDirection::left ? 1 : 2;
    // Don't call this function twice with the same value.
    ASSERT((explored_ & explore_flag) == 0);
    explored_ |= explore_flag;
    Dout(dc::notice, "explored_ is now " << (explored_ == 1 ? "left" : (explored_ == 2 ? "right" : "left|right")));
  }

  bool done() const
  {
    DoutEntering(dc::notice, "LocalExtreme::done() for extreme at " << vertex_sample_ << '.');
    Dout(dc::notice, "explored_ is " << (explored_ == 1 ? "left" : (explored_ == 2 ? "right" : "left|right")));
    return explored_ == 3;
  }
};

} // namespace gradient_descent
