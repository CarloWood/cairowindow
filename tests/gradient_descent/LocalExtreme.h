#pragma once

#include "Sample.h"
#include "Approximation.h"
#include "HorizontalDirection.h"
#include <list>
#include "debug.h"

namespace gradient_descent {

class LocalExtreme
{
 public:
  using extremes_type = std::list<LocalExtreme>;
  static extremes_type s_dummy_list;

 protected:
  Approximation approximation_;         // The (cubic) polynomial approximation around this extreme.
  Sample cp_sample_;                    // Sample taken at the critical point of the approximation_;
  double energy_;                       // The maximum height (Lw) that can be reached from here.
  int explored_{0};                     // Bit 1: exploration to the left of this extreme has started.
                                        // Bit 2: same, on the right.
  std::array<extremes_type::iterator, 2> neighbors_{s_dummy_list.end(), s_dummy_list.end()};
                                        // The extremes on the left (0) and right (1) of this one, if any.

 public:
  LocalExtreme(Sample const& cp_sample, Approximation&& approximation, double energy) :
    // Note cp_sample isn't really moved - we just use std::move to indicate that it is ok to make a copy here.
    cp_sample_(std::move(const_cast<Sample&>(cp_sample))), approximation_(std::move(approximation)), energy_(energy)
  {
    approximation_.set_is_extreme();
  }

  Approximation& approximation()
  {
    return approximation_;
  }

  Approximation const& approximation() const
  {
    return approximation_;
  }

  Sample const& cp_sample() const { return cp_sample_; }

  bool is_minimum() const
  {
    return approximation_.cubic().second_derivative(cp_sample_.w()) > 0.0;
  }

  double energy() const { return energy_; }

  void explored(HorizontalDirection hdirection)
  {
    DoutEntering(dc::notice, "LocalExtreme::explored(" << hdirection << ") for extreme at " << cp_sample_ << '.');
    ASSERT(hdirection != HorizontalDirection::undecided);
    int explore_flag = hdirection == HorizontalDirection::left ? 1 : 2;
    // Don't call this function twice with the same value.
    ASSERT((explored_ & explore_flag) == 0);
    explored_ |= explore_flag;
    Dout(dc::notice, "explored_ is now " << (explored_ == 1 ? "left" : (explored_ == 2 ? "right" : "left|right")));
  }

  bool done() const
  {
    DoutEntering(dc::notice, "LocalExtreme::done() for extreme at " << cp_sample_ << '.');
    Dout(dc::notice, "explored_ is " << (explored_ == 1 ? "left" : (explored_ == 2 ? "right" : "left|right")));
    return explored_ == 3;
  }

  void set_neighbor(HorizontalDirectionToInt direction, extremes_type::iterator neighbor)
  {
#ifdef CWDEBUG
    // Can only assign each element of neighbors_ once, and only with a non-end iterator.
    ASSERT(neighbors_[direction.as_index()] == s_dummy_list.end());
#endif
    neighbors_[direction.as_index()] = neighbor;
  }
  bool has_neighbor(HorizontalDirectionToInt direction) const
  {
    return neighbors_[direction.as_index()] != s_dummy_list.end();
  }
  extremes_type::iterator neighbor(HorizontalDirectionToInt direction) const
  {
    // Call has_neighbor first.
    ASSERT(neighbors_[direction.as_index()] != s_dummy_list.end());
    return neighbors_[direction.as_index()];
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    bool has_left_neighbor = has_neighbor(HorizontalDirection::left);
    bool has_right_neighbor = has_neighbor(HorizontalDirection::right);
    os << "approximation:" << approximation_ <<
        ", cp_sample:" << cp_sample_ <<
        ", energy:" << energy_ <<
        ", explored:" << explored_ <<
        ", neighbors:{";
    if (has_left_neighbor)
      os << "left:" << neighbor(HorizontalDirection::left)->cp_sample().w();
    if (has_right_neighbor)
    {
      if (has_left_neighbor)
        os << ", ";
      os << "right:" << neighbor(HorizontalDirection::right)->cp_sample().w();
    }
    os << '}';
  }
#endif
};

} // namespace gradient_descent
