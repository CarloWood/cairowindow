#pragma once

#include "Weight.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

// This class represents a point for which the value of L(w) and L'(w) have been calculated.
// It is a candidate to be a (local) minimum of L.
//
class Sample : protected Weight
{
 private:
  double Lw_;           // Cache of L(w).
  double dLdw_;         // Cache of L'(w).

 public:
  // Create an uninitialized Sample.
  Sample() = default;

  Sample(Weight w, double Lw, double dLdw) : Weight(w), Lw_(Lw), dLdw_(dLdw) { }

  double w() const
  {
    return w_;
  }

  double Lw() const
  {
    return Lw_;
  }

  double dLdw() const
  {
    return dLdw_;
  }

  void set_values(double w, double Lw, double dLdw)
  {
    w_ = w;
    Lw_ = Lw;
    dLdw_ = dLdw;
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{w:" << w_ << ", Lw:" << Lw_ << ", dLdw:" << dLdw_ << "}";
  }
#endif
};

} // namespace gradient_descent
