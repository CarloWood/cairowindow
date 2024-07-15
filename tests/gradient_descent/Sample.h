#pragma once

#include "Weight.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include <iostream>
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

  // Disallow copying because we use pointers to Sample's all the time.
  // If you know what you're doing, you can use this to make a copy however.
  Sample(Sample&& orig) : Weight(orig), Lw_(orig.Lw_), dLdw_(orig.dLdw_) { }

  // Same for assignment.
  Sample& operator=(Sample&& orig) { w_ = orig.w_; Lw_ = orig.Lw_; dLdw_ = orig.dLdw_; return *this; }

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
