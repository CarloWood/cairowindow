#pragma once

#ifdef CWDEBUG
#include "events/Events.h"
#include "utils/has_print_on.h"

namespace gradient_descent {

using utils::has_print_on::operator<<;

// Event types.

enum event_type
{
  reset_event,
  difference_event,
  fourth_degree_approximation_event,
  derivative_event,
  quotient_event,
  quadratic_polynomial_event,
  kinetic_energy_event
};

class ResetEventData
{
 public:
  ResetEventData() = default;

  void print_on(std::ostream& os) const
  {
    os << "ResetEventData:{}";
  }
};

class DifferenceEventData
{
 protected:
  double w_;
  double expected_Lw_;
  double Lw_;

 public:
  DifferenceEventData(double w, double expected_Lw, double Lw) : w_(w), expected_Lw_(expected_Lw), Lw_(Lw) { }

  double w() const { return w_; }
  double expected_Lw() const { return expected_Lw_; }
  double Lw() const { return Lw_; }

  void print_on(std::ostream& os) const;
};

class PolynomialEventData
{
 private:
  math::Polynomial const& polynomial_;

 public:
  PolynomialEventData(math::Polynomial const& polynomial) : polynomial_(polynomial) { }

  math::Polynomial const& polynomial() const { return polynomial_; }

  void print_on(std::ostream& os) const;
};

class FourthDegreeApproximationEventData : public PolynomialEventData
{
 public:
  using PolynomialEventData::PolynomialEventData;
};

class DerivativeEventData : public PolynomialEventData
{
 public:
  using PolynomialEventData::PolynomialEventData;
};

class QuotientEventData : public PolynomialEventData
{
 public:
  using PolynomialEventData::PolynomialEventData;
};

class QuadraticPolynomialEventData
{
 private:
  math::QuadraticPolynomial const& quadratic_polynomial_;

 public:
  QuadraticPolynomialEventData(math::QuadraticPolynomial const& quadratic_polynomial) : quadratic_polynomial_(quadratic_polynomial) { }

  math::QuadraticPolynomial const& quadratic_polynomial() const { return quadratic_polynomial_; }

  void print_on(std::ostream& os) const;
};

class KineticEnergyEventData
{
 protected:
  double max_Lw_;

 public:
  KineticEnergyEventData(double max_Lw) : max_Lw_(max_Lw) { }

  double max_Lw() const { return max_Lw_; }

  void print_on(std::ostream& os) const
  {
    os << "KineticEnergyEventData:{" << max_Lw_ << "}";
  }
};

class AlgorithmEventData
{
 private:
  std::variant<ResetEventData, DifferenceEventData, FourthDegreeApproximationEventData, QuotientEventData, DerivativeEventData,
    QuadraticPolynomialEventData, KineticEnergyEventData> event_data_;

 public:
  AlgorithmEventData() = default;

  AlgorithmEventData(event_type)
  {
    event_data_.emplace<ResetEventData>();
  }

  AlgorithmEventData(event_type, double w, double expected_Lw, double Lw)
  {
    event_data_.emplace<DifferenceEventData>(w, expected_Lw, Lw);
  }

  AlgorithmEventData(event_type type, math::Polynomial const& polynomial)
  {
    if (type == fourth_degree_approximation_event)
      event_data_.emplace<FourthDegreeApproximationEventData>(polynomial);
    else if (type == quotient_event)
      event_data_.emplace<QuotientEventData>(polynomial);
    else if (type == derivative_event)
      event_data_.emplace<DerivativeEventData>(polynomial);
    else
      ASSERT(false);
  }

  AlgorithmEventData(event_type, math::QuadraticPolynomial const& polynomial)
  {
    event_data_.emplace<QuadraticPolynomialEventData>(polynomial);
  }

  AlgorithmEventData(event_type, double max_Lw)
  {
    event_data_.emplace<KineticEnergyEventData>(max_Lw);
  }

  template<typename T>
  T const& get() const { return std::get<T>(event_data_); }

  template<typename T>
  bool is_a() const { return std::holds_alternative<T>(event_data_); }

  void print_on(std::ostream& os) const;
};

struct AlgorithmEventType : public AlgorithmEventData
{
  using AlgorithmEventData::AlgorithmEventData;
  static constexpr bool one_shot = false;
};

} // namespace gradient_descent
#endif // CWDEBUG
