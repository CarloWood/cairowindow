#pragma once

#include "Expression.h"
#include "precedence.h"

namespace symbolic {

template<Expression E>
class Sin : public ExpressionTag
{
 public:
  using arg_type = E;

  static constexpr precedence s_precedence = precedence::symbol;
  static constexpr auto id_range = E::id_range;

 public:
  static constexpr Sin instance() { return {}; }

#ifdef SYMBOLIC_PRINTING
  static void print_on(std::ostream& os)
  {
    os << "sin(";
    E::print_on(os);
    os << ")";
  }
#endif
};

template<Expression E>
constexpr auto sin(E const&)
{
  return Sin<E>{};
}

template<Expression E>
class Cos : public ExpressionTag
{
 public:
  using arg_type = E;

  static constexpr precedence s_precedence = precedence::symbol;
  static constexpr auto id_range = E::id_range;

 public:
  static constexpr Cos instance() { return {}; }

#ifdef SYMBOLIC_PRINTING
  static void print_on(std::ostream& os)
  {
    os << "cos(";
    E::print_on(os);
    os << ")";
  }
#endif
};

template<Expression E>
constexpr auto cos(E const&)
{
  return Cos<E>{};
}

template<Expression E>
class Log : public ExpressionTag
{
 public:
  using arg_type = E;

  static constexpr precedence s_precedence = precedence::symbol;
  static constexpr auto id_range = E::id_range;

 public:
  static constexpr Log instance() { return {}; }

#ifdef SYMBOLIC_PRINTING
  static void print_on(std::ostream& os)
  {
    os << "log(";
    E::print_on(os);
    os << ")";
  }
#endif
};

template<Expression E>
constexpr auto log(E const&)
{
  return Log<E>{};
}

template<Expression E>
class Atan : public ExpressionTag
{
 public:
  using arg_type = E;

  static constexpr precedence s_precedence = precedence::symbol;
  static constexpr auto id_range = E::id_range;

 public:
  static constexpr Atan instance() { return {}; }

#ifdef SYMBOLIC_PRINTING
  static void print_on(std::ostream& os)
  {
    os << "atan(";
    E::print_on(os);
    os << ")";
  }
#endif
};

template<Expression E>
constexpr auto atan(E const&)
{
  return Atan<E>{};
}

} // namespace symbolic
