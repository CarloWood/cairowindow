#pragma once

#include "Sample.h"
#include "Approximation.h"
#include "debug.h"

namespace gradient_descent {

class LocalExtreme
{
 protected:
  Approximation approximation_;         // The (parabolic) polynomial approximation around this extreme.
  Sample vertex_sample_;                // Sample taken at the vertex of the approximation_;
  double energy_;                       // The maximum height (Lw) that can be reached from here.
  int done_{0};                         // Bit 1: exploration to the left of this extreme has been finished.
                                        // Bit 2: same, on the right.
 public:
  LocalExtreme(Sample const& vertex_sample, Approximation const& approximation, double energy) :
    vertex_sample_(vertex_sample), approximation_(approximation), energy_(energy)
  {
    approximation_.set_is_part_of_extreme();
  }

  Approximation& approximation(HorizontalDirection hdirection)
  {
    DoutEntering(dc::notice, "LocalExtreme::approximation(" << hdirection << ") [" << (void*)this << "]");
    if (hdirection != unknown_horizontal_direction)
    {
      int done_flag = hdirection == left ? 1 : 2;
      // We expect this function to only be called once for each hdirection.
      ASSERT((done_ & done_flag) == 0);
      done_ |= done_flag;
    }
    return approximation_;
  }

  bool is_minimum() const
  {
    return approximation_.parabola()[2] > 0.0;
  }

  Sample const& vertex_sample() const { return vertex_sample_; }

  double energy() const { return energy_; }
  bool done() const { return done_ == 3; }
};

} // namespace gradient_descent
