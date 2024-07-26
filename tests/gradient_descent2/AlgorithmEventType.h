#pragma once

#ifdef CWDEBUG
#include "Scale.h"
#include "../Polynomial.h"
#include "../QuadraticPolynomial.h"
#include "../CubicPolynomial.h"
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
  cubic_polynomial_event,
  kinetic_energy_event,
  scale_draw_event,
  scale_erase_event,
  new_sample_event,
  new_local_extreme_event,
  hdirection_known_event
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

class CubicPolynomialEventData
{
 private:
  math::CubicPolynomial const& cubic_polynomial_;

 public:
  CubicPolynomialEventData(math::CubicPolynomial const& cubic_polynomial) : cubic_polynomial_(cubic_polynomial) { }

  math::CubicPolynomial const& cubic_polynomial() const { return cubic_polynomial_; }

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
  Scale const& scale_;
  math::CubicPolynomial const& old_cubic_;

 public:
  ScaleDrawEventData(ScaleUpdate result, Scale const& scale, math::CubicPolynomial const& old_cubic) :
    result_(result), scale_(scale), old_cubic_(old_cubic) { }

  ScaleUpdate result() const { return result_; }
  Scale const& scale() const { return scale_; }
  math::CubicPolynomial const& old_cubic() const { return old_cubic_; }

  void print_on(std::ostream& os) const
  {
    os << "ScaleDrawEventData:{" << result_ << ", " << scale_ << ", " << old_cubic_ << "}";
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

class NewSampleEventData
{
 protected:
  Sample const& new_sample_;

 public:
  NewSampleEventData(Sample const& new_sample) : new_sample_(new_sample) { }

  Sample const& new_sample() const { return new_sample_; }
  std::string label() const { return std::to_string(new_sample_.label()); }

  void print_on(std::ostream& os) const
  {
    os << "NewSampleEventData:{" << new_sample_ << "}";
  }
};

class NewLocalExtremeEventData
{
 protected:
  Sample const& local_extreme_;
  std::string label_;

 public:
  NewLocalExtremeEventData(Sample const& local_extreme, std::string const& label) :
    local_extreme_(local_extreme), label_(label) { }

  Sample const& local_extreme() const { return local_extreme_; }
  std::string const& label() const { return label_; }

  void print_on(std::ostream& os) const
  {
    os << "NewLocalExtremeEventData:{" << local_extreme_ << ", \"" << label_ << "\"}";
  }
};

class HDirectionKnownEventData
{
 protected:
  Sample const& local_extreme_;
  HorizontalDirection hdirection_;

 public:
  HDirectionKnownEventData(Sample const& local_extreme, HorizontalDirection hdirection) :
    local_extreme_(local_extreme), hdirection_(hdirection) { }

  Sample const& local_extreme() const { return local_extreme_; }
  HorizontalDirection hdirection() const { return hdirection_; }

  void print_on(std::ostream& os) const
  {
    os << "HDirectionKnownEventData:{" << local_extreme_ << ", " << hdirection_ << "}";
  }
};

class AlgorithmEventData
{
 private:
  std::variant<ResetEventData, DifferenceEventData, FourthDegreeApproximationEventData, QuotientEventData, DerivativeEventData,
    QuadraticPolynomialEventData, KineticEnergyEventData, ScaleDrawEventData, ScaleEraseEventData, NewSampleEventData,
    CubicPolynomialEventData, NewLocalExtremeEventData, HDirectionKnownEventData> event_data_;

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

  AlgorithmEventData(event_type type, double w, double expected_Lw, double Lw)
  {
    ASSERT(type == difference_event);
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

  AlgorithmEventData(event_type type, math::QuadraticPolynomial const& polynomial)
  {
    ASSERT(type == quadratic_polynomial_event);
    event_data_.emplace<QuadraticPolynomialEventData>(polynomial);
  }

  AlgorithmEventData(event_type type, math::CubicPolynomial const& polynomial)
  {
    ASSERT(type == cubic_polynomial_event);
    event_data_.emplace<CubicPolynomialEventData>(polynomial);
  }

  AlgorithmEventData(event_type type, double max_Lw)
  {
    ASSERT(type == kinetic_energy_event);
    event_data_.emplace<KineticEnergyEventData>(max_Lw);
  }

  AlgorithmEventData(event_type type, ScaleUpdate result, Scale const& scale, math::CubicPolynomial const& old_cubic)
  {
    ASSERT(type == scale_draw_event);
    event_data_.emplace<ScaleDrawEventData>(result, scale, old_cubic);
  }

  AlgorithmEventData(event_type type, Sample const& current)
  {
    ASSERT(type == new_sample_event);
    event_data_.emplace<NewSampleEventData>(current);
  }

  AlgorithmEventData(event_type type, Sample const& current, std::string const& label)
  {
    ASSERT(type == new_local_extreme_event);
    event_data_.emplace<NewLocalExtremeEventData>(current, label);
  }

  AlgorithmEventData(event_type type, Sample const& current, HorizontalDirection hdirection)
  {
    ASSERT(type == hdirection_known_event);
    event_data_.emplace<HDirectionKnownEventData>(current, hdirection);
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
