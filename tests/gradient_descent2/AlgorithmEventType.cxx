#include "sys.h"
#include "AlgorithmEventType.h"

namespace gradient_descent {

#ifdef CWDEBUG
void DifferenceEventData::print_on(std::ostream& os) const
{
  os << "DifferenceEventData:{w:" << w_ << ", expected_Lw:" << expected_Lw_ << ", Lw:" << Lw_ << "}";
}

void PolynomialEventData::print_on(std::ostream& os) const
{
  os << "PolynomialEventData:{" << polynomial_ << "}";
}

void QuadraticPolynomialEventData::print_on(std::ostream& os) const
{
  os << "QuadraticPolynomialEventData:{" << quadratic_polynomial_ << "}";
}

void CubicPolynomialEventData::print_on(std::ostream& os) const
{
  os << "CubicPolynomialEventData:{" << cubic_polynomial_ << "}";
}

void AlgorithmEventData::print_on(std::ostream& os) const
{
  if (std::holds_alternative<ResetEventData>(event_data_))
    std::get<ResetEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<DifferenceEventData>(event_data_))
    std::get<DifferenceEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<FourthDegreeApproximationEventData>(event_data_))
    std::get<FourthDegreeApproximationEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<DerivativeEventData>(event_data_))
    std::get<DerivativeEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<QuotientEventData>(event_data_))
    std::get<QuotientEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<QuadraticPolynomialEventData>(event_data_))
    std::get<QuadraticPolynomialEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<KineticEnergyEventData>(event_data_))
    std::get<KineticEnergyEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<ScaleDrawEventData>(event_data_))
    std::get<ScaleDrawEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<LeftOfRightOfEventData>(event_data_))
    std::get<LeftOfRightOfEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<JumpPointEventData>(event_data_))
    std::get<JumpPointEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<ScaleEraseEventData>(event_data_))
    std::get<ScaleEraseEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<NewSampleEventData>(event_data_))
    std::get<NewSampleEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<CubicPolynomialEventData>(event_data_))
    std::get<CubicPolynomialEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<NewLocalExtremeEventData>(event_data_))
    std::get<NewLocalExtremeEventData>(event_data_).print_on(os);
  else if (std::holds_alternative<HDirectionKnownEventData>(event_data_))
    std::get<HDirectionKnownEventData>(event_data_).print_on(os);
  else
    // Missing implementation.
    ASSERT(false);
}
#endif

} // namespace gradient_descent
