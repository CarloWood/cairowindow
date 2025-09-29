#pragma once

#ifdef CWDEBUG
#include "HorizontalDirection.h"
#include "ScaleUpdate.h"
#include "SampleNode.h"
#include "math/Polynomial.h"
#include "math/QuadraticPolynomial.h"
#include "math/CubicPolynomial.h"
#include "events/Events.h"
#include "utils/has_print_on.h"
#include <variant>

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
  left_of_right_of_event,
  jump_point_event,
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
  math::Polynomial<double> const& polynomial_;

 public:
  PolynomialEventData(math::Polynomial<double> const& polynomial) : polynomial_(polynomial) { }

  math::Polynomial<double> const& polynomial() const { return polynomial_; }

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
  math::QuadraticPolynomial<double> const& quadratic_polynomial_;

 public:
  QuadraticPolynomialEventData(math::QuadraticPolynomial<double> const& quadratic_polynomial) : quadratic_polynomial_(quadratic_polynomial) { }

  math::QuadraticPolynomial<double> const& quadratic_polynomial() const { return quadratic_polynomial_; }

  void print_on(std::ostream& os) const;
};

class CubicPolynomialEventData
{
 private:
  math::CubicPolynomial<double> const& cubic_polynomial_;
  int index_;

 public:
  CubicPolynomialEventData(math::CubicPolynomial<double> const& cubic_polynomial, HorizontalDirectionToInt index) :
    cubic_polynomial_(cubic_polynomial), index_(index.as_index()) { }

  math::CubicPolynomial<double> const& cubic_polynomial() const { return cubic_polynomial_; }
  int index() const { return index_; }

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
  SampleNode const& sample_node_;
  math::CubicPolynomial<double> const& old_cubic_;

 public:
  ScaleDrawEventData(ScaleUpdate result, SampleNode const& sample_node, math::CubicPolynomial<double> const& old_cubic) :
    result_(result), sample_node_(sample_node), old_cubic_(old_cubic) { }

  ScaleUpdate result() const { return result_; }
  SampleNode const& sample_node() const { return sample_node_; }
  math::CubicPolynomial<double> const& old_cubic() const { return old_cubic_; }

  void print_on(std::ostream& os) const
  {
    os << "ScaleDrawEventData:{" << result_ << ", " << sample_node_ << ", " << old_cubic_ << "}";
  }
};

class LeftOfRightOfEventData
{
 protected:
  SampleNode const* left_of_;
  SampleNode const* right_of_;

 public:
  LeftOfRightOfEventData(SampleNode const* left_of, SampleNode const* right_of) :
    left_of_(left_of), right_of_(right_of) { }

  SampleNode const* left_of() const { return left_of_; }
  SampleNode const* right_of() const { return right_of_; }

  void print_on(std::ostream& os) const
  {
    os << "LeftOfRightOfEventData:{";
    if (right_of_)
      os << "[" << right_of_->label() << "], ";
    else
      os << "-inf, ";
    if (left_of_)
      os << "[" << left_of_->label() << "], ";
    else
      os << "+inf, ";
    os << "}";
  }
};

class JumpPointEventData
{
 protected:
  ExtremeType next_extreme_type_;
  double critical_point_w_;
  double critical_point_Lw_;

 public:
  JumpPointEventData(ExtremeType next_extreme_type, double critical_point_w, double critical_point_Lw) :
    next_extreme_type_(next_extreme_type), critical_point_w_(critical_point_w), critical_point_Lw_(critical_point_Lw) { }

  ExtremeType next_extreme_type() const { return next_extreme_type_; }
  double critical_point_w() const { return critical_point_w_; }
  double critical_point_Lw() const { return critical_point_Lw_; }

  void print_on(std::ostream& os) const
  {
    os << "JumpPointEventData:{";
    os << next_extreme_type_ << ", " << critical_point_w_ << ", " << critical_point_Lw_ << "}";
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
  SampleNode const& extreme_cubic_;
  std::string label_;

 public:
  NewLocalExtremeEventData(SampleNode const& extreme_cubic, std::string const& label) :
    extreme_cubic_(extreme_cubic), label_(label) { }

  SampleNode const& extreme_cubic() const { return extreme_cubic_; }
  std::string const& label() const { return label_; }

  void print_on(std::ostream& os) const
  {
    os << "NewLocalExtremeEventData:{" << extreme_cubic_ << ", \"" << label_ << "\"}";
  }
};

class HDirectionKnownEventData
{
 protected:
  SampleNode const& extreme_cubic_;
  HorizontalDirection hdirection_;

 public:
  HDirectionKnownEventData(SampleNode const& extreme_cubic, HorizontalDirection hdirection) :
    extreme_cubic_(extreme_cubic), hdirection_(hdirection) { }

  SampleNode const& extreme_cubic() const { return extreme_cubic_; }
  HorizontalDirection hdirection() const { return hdirection_; }

  void print_on(std::ostream& os) const
  {
    os << "HDirectionKnownEventData:{" << extreme_cubic_ << ", " << hdirection_ << "}";
  }
};

class AlgorithmEventData
{
 private:
  std::variant<ResetEventData, DifferenceEventData, FourthDegreeApproximationEventData,
    QuotientEventData, DerivativeEventData, QuadraticPolynomialEventData, KineticEnergyEventData, ScaleDrawEventData,
    LeftOfRightOfEventData, JumpPointEventData, ScaleEraseEventData, NewSampleEventData, CubicPolynomialEventData,
    NewLocalExtremeEventData, HDirectionKnownEventData> event_data_;

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

  AlgorithmEventData(event_type type, math::Polynomial<double> const& polynomial)
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

  AlgorithmEventData(event_type type, math::QuadraticPolynomial<double> const& polynomial)
  {
    ASSERT(type == quadratic_polynomial_event);
    event_data_.emplace<QuadraticPolynomialEventData>(polynomial);
  }

  AlgorithmEventData(event_type type, math::CubicPolynomial<double> const& polynomial, HorizontalDirection index)
  {
    ASSERT(type == cubic_polynomial_event);
    event_data_.emplace<CubicPolynomialEventData>(polynomial, index);
  }

  AlgorithmEventData(event_type type, double max_Lw)
  {
    ASSERT(type == kinetic_energy_event);
    event_data_.emplace<KineticEnergyEventData>(max_Lw);
  }

  AlgorithmEventData(event_type type, ScaleUpdate result, SampleNode const& sample_node, math::CubicPolynomial<double> const& old_cubic)
  {
    ASSERT(type == scale_draw_event);
    event_data_.emplace<ScaleDrawEventData>(result, sample_node, old_cubic);
  }

  AlgorithmEventData(event_type type, SampleNode const* left_of, SampleNode const* right_of)
  {
    ASSERT(type == left_of_right_of_event);
    event_data_.emplace<LeftOfRightOfEventData>(left_of, right_of);
  }

  AlgorithmEventData(event_type type, ExtremeType next_extreme_type, double critical_point_w, double critical_point_Lw)
  {
    ASSERT(type == jump_point_event);
    event_data_.emplace<JumpPointEventData>(next_extreme_type, critical_point_w, critical_point_Lw);
  }

  AlgorithmEventData(event_type type, Sample const& current)
  {
    ASSERT(type == new_sample_event);
    event_data_.emplace<NewSampleEventData>(current);
  }

  AlgorithmEventData(event_type type, SampleNode const& current, std::string const& label)
  {
    ASSERT(type == new_local_extreme_event);
    event_data_.emplace<NewLocalExtremeEventData>(current, label);
  }

  AlgorithmEventData(event_type type, SampleNode const& current, HorizontalDirection hdirection)
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
