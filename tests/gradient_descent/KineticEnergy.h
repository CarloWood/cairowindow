#pragma once

#include <cmath>
#include "debug.h"
#ifdef CWDEBUG
#include "AlgorithmEventType.h"
#include "utils/has_print_on.h"
#endif

namespace gradient_descent {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

// A class describing the amount of "kinetic energy" (elsewhere in literature more often known as "momentum")
// that we have, which determines the maximum overshoot uphill that we allow. Having a certain momentum or
// kinetic energy helps in overcoming small bumbs and not get stuck in the first local minimum that we encounter.
class KineticEnergy
{
  static constexpr double friction = 0.0513;    // (1 - friction)^2 = 0.9

 protected:
  double max_Lw_;       // The maximum height that could be reached with the current amount of kinetic energy.
  double Lw_;           // The current height. The kinetic energy is proportional to max_Lw_ - Lw_.

#ifdef CWDEBUG
  events::Server<AlgorithmEventType>& event_server_;
#endif

 private:
  double max_Lw(double new_Lw) const
  {
    // Reduce the maximum height with a fraction of the (vertical) distance traveled.
    return max_Lw_ - friction * std::abs(new_Lw - Lw_);
  }

 public:
  KineticEnergy(double Lw COMMA_CWDEBUG_ONLY(events::Server<AlgorithmEventType>& event_server)) :
    max_Lw_(Lw), Lw_(Lw) COMMA_CWDEBUG_ONLY(event_server_(event_server)) { }

  void set(double max_Lw, double Lw)
  {
    DoutEntering(dc::notice, "KineticEnergy::set(" << max_Lw << ", " << Lw << ")");

    max_Lw_ = max_Lw;
    Lw_ = Lw;
    Dout(dc::notice, "max_Lw_ set to " << max_Lw_ << "; kinetic energy is now " << (max_Lw_ - Lw_));

#ifdef CWDEBUG
    event_server_.trigger(AlgorithmEventType{kinetic_energy_event, max_Lw_});
#endif
  }

  double energy() const
  {
    return max_Lw_;             // This is the total energy.
  }

  // Returns true upon success. If false is returned the update is rejected and should not take place!
  bool maybe_update(double new_Lw)
  {
    DoutEntering(dc::notice, "KineticEnergy::maybe_update(" << new_Lw << ")");

    // Reduce the maximum height with a fraction of this distance.
    double new_max_Lw = max_Lw(new_Lw);

    // Check if this is a request to go higher than the maximum height.
    if (new_Lw > new_max_Lw)
    {
      Dout(dc::notice, "Rejected because max_Lw_ = " << max_Lw_ << ", Lw_ = " << Lw_ << "; new_Lw = " << new_Lw << " > " <<
          "max_Lw_ - friction * std::abs(new_Lw - Lw_) = " << new_max_Lw);
      return false;
    }

    max_Lw_ = new_max_Lw;
    Lw_ = new_Lw;
    Dout(dc::notice, "max_Lw_ set to " << max_Lw_ << "; kinetic energy is now " << (max_Lw_ - Lw_));

#ifdef CWDEBUG
    event_server_.trigger(AlgorithmEventType{kinetic_energy_event, max_Lw_});
#endif
    return true;
  }

  void update(double new_Lw)
  {
    DoutEntering(dc::notice, "KineticEnergy::update(" << new_Lw << ")");

    // Reduce total energy (the maximum height) as usual, but if we need
    // more energy than we have, just reset the total energy to zero.
    max_Lw_ = std::max(new_Lw, max_Lw(new_Lw));
    Lw_ = new_Lw;
    Dout(dc::notice, "max_Lw_ set to " << max_Lw_ << "; kinetic energy is now " << (max_Lw_ - Lw_));

#ifdef CWDEBUG
    event_server_.trigger(AlgorithmEventType{kinetic_energy_event, max_Lw_});
#endif
  }

#if CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{max height:" << max_Lw_ << ", energy:" << (max_Lw_ - Lw_) << "}";
  }
#endif
};

} // namespace gradient_descent
