#include "sys.h"
#include "QuadraticArcLength.h"
#include "utils/derived_from_template.h"

namespace cairowindow::autodiff {

template<typename T>
struct DependentFalse : std::false_type
{
};

template<typename T>
double evaluate(T const&)
{
  static_assert(DependentFalse<T>::value, "Not implemented.");
  AI_NEVER_REACHED
}

template<uint32_t symbol_id_bit, AutoDiffExceptSymbol T>
void differentiate(T const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  static_assert(DependentFalse<T>::value, "Not implemented.");
  AI_NEVER_REACHED
}

//-----------------------------------------------------------------------------
// Zero

template<>
double evaluate(Zero const&)
{
  return 0.0;
}

template<uint32_t symbol_id_bit>
auto differentiate(Zero const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  DoutEntering(dc::notice, "differentiate(" << expression << ", " << symbol << ")");
  return Zero{};
}

//-----------------------------------------------------------------------------
// Constant

template<>
double evaluate(Constant const& arg)
{
  return arg.value();
}

template<uint32_t symbol_id_bit>
auto differentiate(Constant const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  DoutEntering(dc::notice, "differentiate(" << expression << ", " << symbol << ")");
  return Zero{};
}

//-----------------------------------------------------------------------------
// Symbol

template<uint32_t symbol_id_bit>
double evaluate(Symbol<symbol_id_bit> const& arg)
{
  return arg.value();
}

template<uint32_t symbol_id_bit, uint32_t symbol_id_bit2>
auto differentiate(Symbol<symbol_id_bit2> const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  DoutEntering(dc::notice, "differentiate(" << expression << ", " << symbol << ")");
  if constexpr (symbol_id_bit == symbol_id_bit2)
    return Constant{1.0};
  else
    return Zero{};
}

//-----------------------------------------------------------------------------
// DifferentiableSymbol

template<uint32_t symbol_id_bit>
double evaluate(DifferentiableSymbol<symbol_id_bit> const& arg)
{
  return arg.value();
}

//-----------------------------------------------------------------------------

enum Operator
{
  before_plus,
  after_plus,
  before_minus,
  after_minus,
  before_mul,
  after_mul,
  before_div,
  after_div,
  after_negation,
  exponentiated
};

template<Operator op, typename T>
bool needs_parens_(T const& arg)
{
  if constexpr (op == after_negation)
    return true;
  else
    return arg.needs_parens(op);
}

template<Operator op>
bool needs_parens_(Zero const& arg)
{
  return false;
}

template<Operator op>
bool needs_parens_(Constant const& arg)
{
  return false;
}

template<Operator op, uint32_t symbol_id_bit>
bool needs_parens_(Symbol<symbol_id_bit> const& arg)
{
  return false;
}

template<Operator op, uint32_t symbol_id_bit>
bool needs_parens_(DifferentiableSymbol<symbol_id_bit> const& arg)
{
  return false;
}

//-----------------------------------------------------------------------------
// Negation

template<AutoDiff T1>
class Negation : public AutoDiffExpressionTag
{
 private:
  T1 value_;

 public:
  Negation(T1 const& value) : value_(value) { }

  T1 const& value() const { return value_; }

  T1 const& operator-() const { return value_; }

#ifdef CWDEBUG
  bool needs_parens(Operator op) const
  {
    return false;
  }

  void print_on(std::ostream& os) const
  {
    bool need_parens = needs_parens_<after_negation>(value_);
    os << "-!";
    if (need_parens)
      os << '(';
    value_.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<typename T>
concept AutoDiffExceptNegation = std::is_base_of_v<AutoDiffTag, T> && !utils::derived_from_template_v<T, Negation>;

template<typename T>
concept AutoDiffExceptNegationOrZero = AutoDiffExceptZero<T> && !utils::derived_from_template_v<T, Negation>;

template<AutoDiff T1>
double evaluate(Negation<T1> const& negation)
{
  return -evaluate(negation.value());
}

template<uint32_t symbol_id_bit>
Negation<Symbol<symbol_id_bit>> operator-(Symbol<symbol_id_bit> const& arg)
{
  return arg;
}

template<uint32_t symbol_id_bit>
Negation<DifferentiableSymbol<symbol_id_bit>> operator-(DifferentiableSymbol<symbol_id_bit> const& arg)
{
  return arg;
}

template<uint32_t symbol_id_bit, AutoDiff T>
auto differentiate(Negation<T> const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  DoutEntering(dc::notice, "differentiate(" << expression << ", " << symbol << ")");
  return -differentiate(expression.value(), symbol);
}

//-----------------------------------------------------------------------------
// Addition

template<AutoDiff T1, AutoDiff T2>
class Addition : public AutoDiffExpressionTag
{
 private:
  T1 arg1_;
  T2 arg2_;

 public:
  Addition(T1 const& arg1, T2 const& arg2) : arg1_(arg1), arg2_(arg2) { }

  T1 const& arg1() const { return arg1_; }
  T2 const& arg2() const { return arg2_; }

#ifdef CWDEBUG
  bool needs_parens(Operator op) const
  {
    switch (op)
    {
      // x - (y + z); and all mul/div.
      case before_plus:
      case after_plus:
      case before_minus:
        return false;
      default:
        break;
    }
    return true;
  }

  void print_on(std::ostream& os) const
  {
    bool need_parens = needs_parens_<before_plus>(arg1_);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " + ";
    need_parens = needs_parens_<after_plus>(arg2_);
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<AutoDiff T1, AutoDiff T2>
double evaluate(Addition<T1, T2> const& addition)
{
  return evaluate(addition.arg1()) + evaluate(addition.arg2());
}

template<AutoDiffExceptNegation T1, AutoDiffExceptNegation T2>
Addition<T1, T2> operator+(T1 const& arg1, T2 const& arg2)
{
  return {arg1, arg2};
}

template<AutoDiff T1>
T1 operator+(T1 const& arg1, Zero const&)
{
  return arg1;
}

template<AutoDiffExceptZero T1>
T1 operator+(Zero const&, T1 const& arg2)
{
  return arg2;
}

Constant operator+(Constant const& arg1, Constant const& arg2)
{
  return {arg1.value() + arg2.value()};
}

template<uint32_t symbol_id_bit, AutoDiff T1, AutoDiff T2>
auto differentiate(Addition<T1, T2> const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  DoutEntering(dc::notice, "differentiate(" << expression << ", " << symbol << ")");
  return differentiate(expression.arg1(), symbol) + differentiate(expression.arg2(), symbol);
}

//-----------------------------------------------------------------------------
// Subtraction

template<AutoDiff T1, AutoDiff T2>
class Subtraction : public AutoDiffExpressionTag
{
 private:
  T1 arg1_;
  T2 arg2_;

 public:
  Subtraction(T1 const& arg1, T2 const& arg2) : arg1_(arg1), arg2_(arg2) { }

  T1 const& arg1() const { return arg1_; }
  T2 const& arg2() const { return arg2_; }

#ifdef CWDEBUG
  bool needs_parens(Operator op) const
  {
    switch (op)
    {
      // x - (y - z); and all mul/div.
      case before_plus:
      case after_plus:
      case before_minus:
        return false;
      default:
        break;
    }
    return true;
  }

  void print_on(std::ostream& os) const
  {
    bool need_parens = needs_parens_<before_minus>(arg1_);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " - ";
    need_parens = needs_parens_<after_minus>(arg2_);
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<AutoDiff T1, AutoDiff T2>
double evaluate(Subtraction<T1, T2> const& addition)
{
  return evaluate(addition.arg1()) - evaluate(addition.arg2());
}

template<AutoDiffExceptNegationOrZero T1, AutoDiff T2>
Subtraction<T1, T2> operator+(T1 const& arg1, Negation<T2> const& arg2)
{
  return {arg1, -arg2};
}

template<AutoDiff T1, AutoDiffExceptNegationOrZero T2>
Subtraction<T2, T1> operator+(Negation<T1> const& arg1, T2 const& arg2)
{
  return {arg2, -arg1};
}

template<AutoDiff T1, AutoDiff T2>
auto operator+(Negation<T1> const& arg1, Negation<T2> const& arg2)
{
  return Subtraction{arg1, -arg2};
}

template<AutoDiff T1, AutoDiffExceptNegation T2>
Subtraction<T1, T2> operator-(T1 const& arg1, T2 const& arg2)
{
  return {arg1, arg2};
}

template<AutoDiff T1>
T1 operator-(T1 const& arg1, Zero const&)
{
  return arg1;
}

template<AutoDiffExceptNegationOrZero T1>
Negation<T1> operator-(Zero const&, T1 const& arg2)
{
  return arg2;
}

template<AutoDiff T1>
T1 operator-(Zero const&, Negation<T1> const& arg2)
{
  return -arg2;
}

Constant operator-(Zero const&, Constant const& arg2)
{
  return -arg2.value();
}

Constant operator-(Constant const& arg1, Constant const& arg2)
{
  return arg1.value() - arg2.value();
}

template<AutoDiffExceptNegation T1, AutoDiff T2>
Addition<T1, T2> operator-(T1 const& arg1, Negation<T2> const& arg2)
{
  return {arg1, -arg2};
}

template<AutoDiff T1>
Subtraction<Constant, T1> operator-(Negation<T1> const& arg1, Constant const& arg2)
{
  return {-arg2, -arg1};
}

template<AutoDiff T1, AutoDiff T2>
Subtraction<T2, T1> operator-(Negation<T1> const& arg1, Negation<T2> const& arg2)
{
  return {-arg2, -arg1};
}

template<uint32_t symbol_id_bit, AutoDiff T1, AutoDiff T2>
auto differentiate(Subtraction<T1, T2> const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  DoutEntering(dc::notice, "differentiate(" << expression << ", " << symbol << ")");
  return differentiate(expression.arg1(), symbol) - differentiate(expression.arg2(), symbol);
}

//-----------------------------------------------------------------------------
// Multiplication

template<AutoDiff T1, AutoDiff T2>
requires (!std::is_same_v<T2, Constant>)
class Multiplication : public AutoDiffExpressionTag
{
 public:
  using arg1_type = T1;
  using arg2_type = T2;

 private:
  T1 arg1_;
  T2 arg2_;

 public:
  Multiplication(T1 const& arg1, T2 const& arg2) : arg1_(arg1), arg2_(arg2) { }

  T1 const& arg1() const { return arg1_; }
  T2 const& arg2() const { return arg2_; }

#ifdef CWDEBUG
  bool needs_parens(Operator op) const
  {
    //FIXME: remove this
    return true;
    // x / (y * z)
//    return op == after_div || op == exponentiated;
  }

  void print_on(std::ostream& os) const
  {
    bool need_parens = needs_parens_<before_mul>(arg1_);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " * ";
    need_parens = needs_parens_<after_mul>(arg2_);
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<AutoDiff T1, AutoDiff T2>
double evaluate(Multiplication<T1, T2> const& addition)
{
  return evaluate(addition.arg1()) * evaluate(addition.arg2());
}

template<AutoDiff T1, AutoDiff T2>
Multiplication<T1, T2> operator*(T1 const& arg1, T2 const& arg2)
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return {arg1, arg2};
}

// Let A and B be constants (where C = A * B), then we want to reorder multiplications as follows:
// M(A, M(B, e)) --> M(M(A, B), e) --> M(C, e)
template<AutoDiff T3>
auto operator*(Constant const& arg1, Multiplication<Constant, T3> const& arg2)
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return Multiplication{arg1 * arg2.arg1(), arg2.arg2()};
}

// M(M(A, e1), M(B, e2)) --> M(M(A, B), M(e1, e2)) --> M(C, M(e1, e2))
template<AutoDiff T2, AutoDiff T4>
auto operator*(Multiplication<Constant, T2> const& arg1, Multiplication<Constant, T4> const& arg2)
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return Multiplication{arg1.arg1() * arg2.arg1(), arg1.arg2() * arg2.arg2()};
}

// M(M(A, e1), e2) --> M(A, M(e1, e2))
template<AutoDiff T2, AutoDiffExceptZero T3>
requires (!std::is_same_v<T3, Constant>)
auto operator*(Multiplication<Constant, T2> const& arg1, T3 const& arg2)
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return Multiplication{arg1.arg1(), arg1.arg2() * arg2};
}

template<AutoDiffExceptZero T1>
auto operator*(T1 const& arg1, Constant const& arg2)
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return arg2 * arg1;
}

template<AutoDiff T1>
Zero operator*(T1 const& CWDEBUG_ONLY(arg1), Zero const& CWDEBUG_ONLY(arg2))
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return {};
}

template<AutoDiffExceptZero T1>
Zero operator*(Zero const& CWDEBUG_ONLY(arg1), T1 const& CWDEBUG_ONLY(arg2))
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return {};
}

Constant operator*(Constant const& arg1, Constant const& arg2)
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return arg1.value() * arg2.value();
}

template<AutoDiff T2>
Multiplication<Constant, T2> operator*(Constant const& arg1, Negation<T2> const& arg2)
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return {-arg1, -arg2};
}

template<AutoDiff T1>
Multiplication<Constant, T1> operator*(Negation<T1> const& arg1, Constant const& arg2)
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return {-arg2, -arg1};
}

template<AutoDiff T1, AutoDiff T2>
auto operator*(Negation<T1> const& arg1, Negation<T2> const& arg2)
{
  DoutEntering(dc::notice, "operator*(" << arg1 << ", " << arg2 << ")");
  return Multiplication{-arg1, -arg2};
}

template<uint32_t symbol_id_bit, AutoDiff T1, AutoDiff T2>
auto differentiate(Multiplication<T1, T2> const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  DoutEntering(dc::notice, "differentiate(" << expression << ", " << symbol << ")");
  return differentiate(expression.arg1(), symbol) * expression.arg2() + expression.arg1() * differentiate(expression.arg2(), symbol);
}

//-----------------------------------------------------------------------------
// Division

template<AutoDiff T1, AutoDiff T2>
class Division : public AutoDiffExpressionTag
{
 private:
  T1 arg1_;
  T2 arg2_;

 public:
  Division(T1 const& arg1, T2 const& arg2) : arg1_(arg1), arg2_(arg2) { }

  T1 const& arg1() const { return arg1_; }
  T2 const& arg2() const { return arg2_; }

#ifdef CWDEBUG
  bool needs_parens(Operator op) const
  {
    // x / (y / z)
    return op == after_div || op == exponentiated;
  }

  void print_on(std::ostream& os) const
  {
    bool need_parens = needs_parens_<before_div>(arg1_);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " / ";
    need_parens = needs_parens_<after_div>(arg2_);
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<AutoDiff T1, AutoDiff T2>
double evaluate(Division<T1, T2> const& addition)
{
  return evaluate(addition.arg1()) / evaluate(addition.arg2());
}

template<AutoDiff T1, AutoDiffExceptNegationOrZero T2>
Division<T1, T2> operator/(T1 const& arg1, T2 const& arg2)
{
  DoutEntering(dc::notice, "operator/(" << arg1 << ", " << arg2 << ")");
  return {arg1, arg2};
}

template<AutoDiffExceptZero T1>
Zero operator/(Zero const& CWDEBUG_ONLY(arg1), T1 const& CWDEBUG_ONLY(arg2))
{
  DoutEntering(dc::notice, "operator/(" << arg1 << ", " << arg2 << ")");
  return {};
}

Constant operator/(Constant const& arg1, Constant const& arg2)
{
  DoutEntering(dc::notice, "operator/(" << arg1 << ", " << arg2 << ")");
  return arg1.value() / arg2.value();
}

template<AutoDiffExceptZero T1, AutoDiff T2>
auto operator/(T1 const& arg1, Negation<T2> const& arg2)
{
  DoutEntering(dc::notice, "operator/(" << arg1 << ", " << arg2 << ")");
  return Division{-arg1, -arg2};
}

template<AutoDiff T1>
auto operator/(Negation<T1> const& arg1, Constant const& arg2)
{
  DoutEntering(dc::notice, "operator/(" << arg1 << ", " << arg2 << ")");
  return Division{-arg1, -arg2};
}

//-----------------------------------------------------------------------------
// Exponentiation

template<AutoDiff T1>
class Exponentiation : public AutoDiffExpressionTag
{
 private:
  T1 arg1_;
  int exponent_;

 public:
  Exponentiation(T1 const& arg1, int exponent) : arg1_(arg1), exponent_(exponent)
  {
    ASSERT(exponent < 0 || exponent > 1);
  }

  T1 const& arg1() const { return arg1_; }
  int exponent() const { return exponent_; }

#ifdef CWDEBUG
  bool needs_parens(Operator op) const
  {
    return op == exponentiated;
  }

  void print_on(std::ostream& os) const
  {
    bool need_parens = needs_parens_<exponentiated>(arg1_);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    int e = exponent_;
    if (exponent_ < 0)
    {
      os << "⁻";
      e = -exponent_;
    }
    static std::array<char const*, 10> superscript_digit = { "⁰", "¹", "²", "³", "⁴", "⁵", "⁶", "⁷", "⁸", "⁹" };
    int div = 1;
    int e2 = e;
    while (e2 >= 10)
    {
      e2 /= 10;
      div *= 10;
    }
    while (div > 0)
    {
      os << superscript_digit[e / div];
      e %= div;
      div /= 10;
    }
  }
#endif
};

// Division
template<uint32_t symbol_id_bit, AutoDiff T1, AutoDiff T2>
auto differentiate(Division<T1, T2> const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  DoutEntering(dc::notice, "differentiate(" << expression << ", " << symbol << ")");
  return (differentiate(expression.arg1(), symbol) * expression.arg2() -
      expression.arg1() * differentiate(expression.arg2(), symbol)) / Exponentiation{expression.arg2(), 2};
}

template<AutoDiff T1>
double evaluate(Exponentiation<T1> const& exponentiation)
{
  return std::pow(evaluate(exponentiation.arg1()), exponentiation.exponent());
}

template<AutoDiff T1>
auto operator^(T1 const& arg, int exponent)
{
  DoutEntering(dc::notice, "operator^(" << arg << ", " << exponent << ")");
  return Exponentiation{arg, exponent};
}

template<AutoDiff T1>
auto operator^(Exponentiation<T1> const& arg, int exponent)
{
  DoutEntering(dc::notice, "operator^(" << arg << ", " << exponent << ")");
  return Exponentiation{arg.arg1(), arg.exponent() * exponent};
}

template<uint32_t symbol_id_bit, AutoDiff T1>
auto differentiate(Exponentiation<T1> const& expression, DifferentiableSymbol<symbol_id_bit> const& symbol)
{
  DoutEntering(dc::notice, "differentiate(" << expression << ", " << symbol << ")");
  return Constant{static_cast<double>(expression.exponent())} * Exponentiation{expression.arg1(), expression.exponent() - 1} *
    differentiate(expression.arg1(), symbol);
}

//-----------------------------------------------------------------------------

// After calling BezierCurve::quadratic_from(double v0qa, double v1qa) (and that returned true),
// the following functions can be called to (redo) the calculations for the arc length;

// V₀_x.
double QuadraticArcLength::v0x(double v0qa, double v1qa)
{
  // Q1 and N1 are constants (do not depend on v0qa or v1qa).
  double Q1x = m_.coefficient[2].x() + m_.coefficient[1].x();
  double N1x = -(m_.coefficient[2].y() + m_.coefficient[1].y());
  double v0_div_q1 = 2.0 * std::sin(v1qa) / std::sin(v1qa - v0qa);
  return v0_div_q1 * (std::cos(v0qa) * Q1x + std::sin(v0qa) * N1x);
}

// V₀_y.
double QuadraticArcLength::v0y(double v0qa, double v1qa)
{
  double Q1y = m_.coefficient[2].y() + m_.coefficient[1].y();
  double N1y = m_.coefficient[2].x() + m_.coefficient[1].x();
  double v0_div_q1 = 2.0 * std::sin(v1qa) / std::sin(v1qa - v0qa);
  return v0_div_q1 * (std::cos(v0qa) * Q1y + std::sin(v0qa) * N1y);
}

// |V₀|² = V₀_x² + V₀_y².
double QuadraticArcLength::v02()
{
  return V0().length_squared();
}

// A₀_x = 2 (Q1_x - V₀_x).
double QuadraticArcLength::A0x()
{
  double Q1x = m_.coefficient[2].x() + m_.coefficient[1].x();
  return 2.0 * (Q1x - V0().x());
}

// A₀_y = 2 (Q1_y - V₀_y).
double QuadraticArcLength::A0y()
{
  double Q1y = m_.coefficient[2].y() + m_.coefficient[1].y();
  return 2.0 * (Q1y - V0().y());
}

// |A₀|² = A₀_x² + A₀_y².
double QuadraticArcLength::a02()
{
  return A0().length_squared();
}

// v0 = |V₀| = √(V₀_x² + V₀_y²).
double QuadraticArcLength::v0()
{
  return std::sqrt(v02());
}

// a0 = |A₀| = √(A₀_x² + A₀_y²) = 2√((Q1_x - V₀_x)² + (Q1_y - V₀_y)²).
double QuadraticArcLength::a0()
{
  return std::sqrt(a02());
}

// z = V₀·A₀ = V₀_x A₀_x + V₀_y A₀_y = V₀_x (2 (Q1_x - V₀_x)) + V₀_y (2 (Q1_y - V₀_y)) =
//     2 (V₀_x Q1_x - V₀_x² + V₀_y Q1_y - V₀_y²).
double QuadraticArcLength::z()
{
  return V0().dot(A0());
}

double QuadraticArcLength::s()
{
  return std::sqrt(v02() + 2.0 * z() + a02());
}

// a03 = |A₀|³
double QuadraticArcLength::a03()
{
  return a02() * a0();
}

// (z * a0 + a03) * s
double QuadraticArcLength::za0pa03s()
{
  return (z() * a0() + a03()) * s();
}

double QuadraticArcLength::za0v0()
{
  return z() * a0() * v0();
}

double QuadraticArcLength::v02a02mz2()
{
  return v02() * a02() - utils::square(z());
}

double QuadraticArcLength::zpa02pa0s()
{
  return z() + a02() + a0() * s();
}

double QuadraticArcLength::zpv0a0()
{
  return z() + v0() * a0();
}

double QuadraticArcLength::the_log()
{
  return std::log(zpa02pa0s() / zpv0a0());
}

double QuadraticArcLength::the_enumerator()
{
  return za0pa03s() - za0v0() + v02a02mz2() * the_log();
}

double QuadraticArcLength::q_arc_length()
{
  return the_enumerator() / (2.0 * a03());
}

// Test symbols.
namespace symbol {
static constexpr uint32_t x = 100;
static constexpr uint32_t y = 101;
} // namespace symbol

void QuadraticArcLength::test()
{
  Zero zero1;
  Zero zero2;

  Dout(dc::notice, "zero = " << zero1);
  Dout(dc::notice, "evaluate(zero) = " << evaluate(zero1));
  Dout(dc::notice, "-zero = " << (-zero1) << " = " << evaluate(-zero1));

  Dout(dc::notice, "zero + zero = " << (zero1 + zero2));
  Dout(dc::notice, "zero - zero = " << (zero1 - zero2));
  Dout(dc::notice, "zero * zero = " << (zero1 * zero2));

  Constant constant1(42);

  Dout(dc::notice, "constant = " << constant1);
  Dout(dc::notice, "evaluate(constant1) = " << evaluate(constant1));
  Dout(dc::notice, "-constant = " << (-constant1) << " = " << evaluate(-constant1));

  Dout(dc::notice, "zero + constant(" << evaluate(constant1) << ") = " << (zero1 + constant1));
  Dout(dc::notice, "zero - constant(" << evaluate(constant1) << ") = " << (zero1 - constant1));
  Dout(dc::notice, "zero * constant(" << evaluate(constant1) << ") = " << (zero1 * constant1));
  Dout(dc::notice, "zero / constant(" << evaluate(constant1) << ") = " << (zero1 / constant1));

  Dout(dc::notice, "constant(" << evaluate(constant1) << ") + zero = " << (constant1 + zero1));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") - zero = " << (constant1 - zero1));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") * zero = " << (constant1 * zero1));

  Constant constant2(13);

  Dout(dc::notice, "constant(" << evaluate(constant1) << ") + constant(" << evaluate(constant2) << ") = " << (constant1 + constant2));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") - constant(" << evaluate(constant2) << ") = " << (constant1 - constant2));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") * constant(" << evaluate(constant2) << ") = " << (constant1 * constant2));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") / constant(" << evaluate(constant2) << ") = " << (constant1 / constant2));

  DifferentiableSymbol<symbol::x> x("x");
  x = 100;

  Dout(dc::notice, "x = " << x);
  Dout(dc::notice, "evaluate(x) = " << evaluate(x));

  Dout(dc::notice, "zero + x = " << (zero1 + x) << " = " << evaluate(zero1 + x));
  Dout(dc::notice, "zero - x = " << (zero1 - x) << " = " << evaluate(zero1 - x));
  Dout(dc::notice, "zero * x = " << (zero1 * x) << " = " << evaluate(zero1 * x));
  Dout(dc::notice, "zero / x = " << (zero1 / x) << " = " << evaluate(zero1 / x));

  Dout(dc::notice, "x + zero = " << (x + zero1) << " = " << evaluate(x + zero1));
  Dout(dc::notice, "x - zero = " << (x - zero1) << " = " << evaluate(x - zero1));
  Dout(dc::notice, "x * zero = " << (x * zero1) << " = " << evaluate(x * zero1));

  Dout(dc::notice, "constant(" << evaluate(constant1) << ") + x = " << (constant1 + x) << " = " << evaluate(constant1 + x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") - x = " << (constant1 - x) << " = " << evaluate(constant1 - x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") * x = " << (constant1 * x) << " = " << evaluate(constant1 * x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") / x = " << (constant1 / x) << " = " << evaluate(constant1 / x));

  Dout(dc::notice, "x + constant(" << evaluate(constant1) << ") = " << (x + constant1) << " = " << evaluate(x + constant1));
  Dout(dc::notice, "x - constant(" << evaluate(constant1) << ") = " << (x - constant1) << " = " << evaluate(x - constant1));
  Dout(dc::notice, "x * constant(" << evaluate(constant1) << ") = " << (x * constant1) << " = " << evaluate(x * constant1));
  Dout(dc::notice, "x / constant(" << evaluate(constant1) << ") = " << (x / constant1) << " = " << evaluate(x / constant1));

  Dout(dc::notice, "-x = " << (-x) << " = " << evaluate(-x));

  Dout(dc::notice, "zero + -x = " << (zero1 + -x) << " = " << evaluate(zero1 + -x));
  Dout(dc::notice, "zero - -x = " << (zero1 - -x) << " = " << evaluate(zero1 - -x));
  Dout(dc::notice, "zero * -x = " << (zero1 * -x) << " = " << evaluate(zero1 * -x));
  Dout(dc::notice, "zero / -x = " << (zero1 / -x) << " = " << evaluate(zero1 / -x));

  Dout(dc::notice, "-x + zero = " << (-x + zero1) << " = " << evaluate(-x + zero1));
  Dout(dc::notice, "-x - zero = " << (-x - zero1) << " = " << evaluate(-x - zero1));
  Dout(dc::notice, "-x * zero = " << (-x * zero1) << " = " << evaluate(-x * zero1));

  Dout(dc::notice, "constant(" << evaluate(constant1) << ") + -x = " << (constant1 + -x) << " = " << evaluate(constant1 + -x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") - -x = " << (constant1 - -x) << " = " << evaluate(constant1 - -x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") * -x = " << (constant1 * -x) << " = " << evaluate(constant1 * -x));
  Dout(dc::notice, "constant(" << evaluate(constant1) << ") / -x = " << (constant1 / -x) << " = " << evaluate(constant1 / -x));

  Dout(dc::notice, "-x + constant(" << evaluate(constant1) << ") = " << (-x + constant1) << " = " << evaluate(-x + constant1));
  Dout(dc::notice, "-x - constant(" << evaluate(constant1) << ") = " << (-x - constant1) << " = " << evaluate(-x - constant1));
  Dout(dc::notice, "-x * constant(" << evaluate(constant1) << ") = " << (-x * constant1) << " = " << evaluate(-x * constant1));
  Dout(dc::notice, "-x / constant(" << evaluate(constant1) << ") = " << (-x / constant1) << " = " << evaluate(-x / constant1));

  Symbol<symbol::y> y("y");
  y = 200;

  Dout(dc::notice, "x + y = " << (x + y) << " = " << evaluate(x + y));
  Dout(dc::notice, "x - y = " << (x - y) << " = " << evaluate(x - y));
  Dout(dc::notice, "x * y = " << (x * y) << " = " << evaluate(x * y));
  Dout(dc::notice, "x / y = " << (x / y) << " = " << evaluate(x / y));

  Dout(dc::notice, "-x + y = " << (-x + y) << " = " << evaluate(-x + y));
  Dout(dc::notice, "-x - y = " << (-x - y) << " = " << evaluate(-x - y));
  Dout(dc::notice, "-x * y = " << (-x * y) << " = " << evaluate(-x * y));
  Dout(dc::notice, "-x / y = " << (-x / y) << " = " << evaluate(-x / y));

  Dout(dc::notice, "x + -y = " << (x + -y) << " = " << evaluate(x + -y));
  Dout(dc::notice, "x - -y = " << (x - -y) << " = " << evaluate(x - -y));
  Dout(dc::notice, "x * -y = " << (x * -y) << " = " << evaluate(x * -y));
  Dout(dc::notice, "x / -y = " << (x / -y) << " = " << evaluate(x / -y));

  Dout(dc::notice, "-x + -y = " << (-x + -y) << " = " << evaluate(-x + -y));
  Dout(dc::notice, "-x - -y = " << (-x - -y) << " = " << evaluate(-x - -y));
  Dout(dc::notice, "-x * -y = " << (-x * -y) << " = " << evaluate(-x * -y));
  Dout(dc::notice, "-x / -y = " << (-x / -y) << " = " << evaluate(-x / -y));

  Dout(dc::notice, differentiate(constant1, x));

  Dout(dc::notice, "∂x/∂x = " << differentiate(x, x));
  Dout(dc::notice, "∂y/∂x = " << differentiate(y, x));

  Dout(dc::notice, "∂-x/∂x = " << differentiate(-x, x));
  Dout(dc::notice, "∂(x + y + 42)/∂x = " << differentiate(x + y + constant1, x));
  Dout(dc::notice, "∂(x + y - 42)/∂x = " << differentiate(x + y - constant1, x));
  Dout(dc::notice, "∂(y - x)/∂x = " << differentiate(y - x, x));

  Dout(dc::notice, "∂(x * x * x)/∂x = " << differentiate(x * x * x, x));

  Dout(dc::notice, "∂((y - x^2)^5)/∂x = " << differentiate((y - (x^2))^5, x));
}

} // namespace cairowindow::autodiff
