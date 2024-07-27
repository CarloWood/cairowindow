#pragma once

#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include "utils/UniqueID.h"
#include <iostream>
#endif
#include "debug.h"

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

// Stores an actual "measurement" of the function under test: L(w).
// Note that this function is defined by the user: we only get samples
// by means of calls to Algorithm::operator().
class Sample
{
 private:
  double w_;            // w value.
  double Lw_;           // The corresponding function value L(w).
  double dLdw_;         // The corresponding derivative value dL/dw at w.
#ifdef CWDEBUG
  utils::UniqueID<int> label_;
#endif

 public:
  Sample(double w, double Lw, double dLdw COMMA_CWDEBUG_ONLY(utils::UniqueIDContext<int>& label_context)) :
    w_(w), Lw_(Lw), dLdw_(dLdw) COMMA_CWDEBUG_ONLY(label_(label_context.get_id())) { }

  // Accessors.
  double w() const { return w_; }
  double Lw() const { return Lw_; }
  double dLdw() const { return dLdw_; }
#ifdef CWDEBUG
  int label() const { return label_; }
#endif

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{w:" << w_ << ", Lw:" << Lw_ << ", dLdw:" << dLdw_ << "}";
  }
#endif
};

} // namespace gradient_descent
