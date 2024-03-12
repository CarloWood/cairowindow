#pragma once

#include "IdRange.h"
#include "Expression.h"
#include "counter.h"
#include <cstring>
#ifdef SYMBOLIC_PRINTING
#include "utils/has_print_on.h"
#include <iostream>
#endif

namespace symbolic {
#ifdef SYMBOLIC_PRINTING
using utils::has_print_on::operator<<;
#endif

template<int Id>
class Symbol : public ExpressionTag
{
 public:
  static constexpr precedence s_precedence = precedence::symbol;
  static constexpr IdRange<Id, Id + 1> id_range;

 private:
  char const* const name_;
  int value_;

 public:
  constexpr Symbol(char const* name) : name_(name), value_(0.0) { }

  constexpr int id() const { return Id; }
  char const* name() const { return name_; }

  Symbol const& operator=(int value)
  {
    value_ = value;
    return *this;
  }

#ifdef SYMBOLIC_PRINTING
  constexpr bool needs_parens(before_or_after position, precedence prec) const { return prec < s_precedence; }
  constexpr bool needs_parens(precedence prec) const { return prec < s_precedence; }

  void print_on(std::ostream& os) const
  {
    os << name_;
  }
#endif
};

template<typename T>
struct is_symbol : std::false_type { };

template<typename T>
constexpr bool is_symbol_v = is_symbol<T>::value;

template<int Id>
struct is_symbol<Symbol<Id>> : std::true_type { };

template<typename E>
concept SymbolType = is_symbol_v<E>;

template<auto Id = int{}, typename T = decltype([]{})>
constexpr auto make_symbol(char const* name)
{
  return Symbol<metahack::unique_id<Id, T>()>(name);
}

template<int Id1, int Id2>
constexpr bool is_same_expression(Symbol<Id1> const& arg1, Symbol<Id2> const& arg2) { return arg1.id() == arg2.id(); }

} // namespace symbolic
