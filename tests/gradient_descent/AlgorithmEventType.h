#pragma once

#ifdef CWDEBUG
#include "HistoryIndex.h"
#include "Scale.h"
#include "../Polynomial.h"
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
  kinetic_energy_event,
  scale_draw_event,
  scale_erase_event,
  history_add_event
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

class ScaleDrawEventData
{
 protected:
  ScaleUpdate result_;
  double x1_;
  double x2_;
  math::QuadraticPolynomial const& old_parabola_;

 public:
  ScaleDrawEventData(ScaleUpdate result, double x1, double x2, math::QuadraticPolynomial const& old_parabola) :
    result_(result), x1_(x1), x2_(x2), old_parabola_(old_parabola) { }

  ScaleUpdate result() const { return result_; }
  double x1() const { return x1_; }
  double x2() const { return x2_; }
  math::QuadraticPolynomial const& old_parabola() const { return old_parabola_; }

  void print_on(std::ostream& os) const
  {
    os << "ScaleDrawEventData:{" << result_ << ", " << x1_ << ", " << x2_ << ", " << old_parabola_ << "}";
  }
};

class ScaleEraseEventData
{
 public:
  void print_on(std::ostream& os) const
  {
    os << "ScaleEraseEventData:{}";
  }
};

class HistoryAddEventData
{
 protected:
  HistoryIndex index_;
  Sample const& current_;
  std::string label_;

 public:
  HistoryAddEventData(HistoryIndex index, Sample const& current, std::string const& label) :
    index_(index), current_(current), label_(label) { }

  HistoryIndex index() const { return index_; }
  Sample const& current() const { return current_; }
  std::string const& label() const { return label_; }

  void print_on(std::ostream& os) const
  {
    os << "HistoryAddEventData:{" << index_ << ", " << current_ << ", \"" << label_ << "\"}";
  }
};

class AlgorithmEventData
{
 private:
  std::variant<ResetEventData, DifferenceEventData, FourthDegreeApproximationEventData, QuotientEventData, DerivativeEventData,
    QuadraticPolynomialEventData, KineticEnergyEventData, ScaleDrawEventData, ScaleEraseEventData, HistoryAddEventData> event_data_;

 public:
  AlgorithmEventData() = default;

  AlgorithmEventData(event_type type)
  {
    if (type == reset_event)
      event_data_.emplace<ResetEventData>();
    else if (type == scale_erase_event)
      event_data_.emplace<ScaleEraseEventData>();
    else
      ASSERT(false);
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

  AlgorithmEventData(event_type, ScaleUpdate result, double x1, double x2, math::QuadraticPolynomial const& old_parabola)
  {
    event_data_.emplace<ScaleDrawEventData>(result, x1, x2, old_parabola);
  }

  AlgorithmEventData(event_type, HistoryIndex index, Sample const& current, std::string const& label)
  {
    event_data_.emplace<HistoryAddEventData>(index, current, label);
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
