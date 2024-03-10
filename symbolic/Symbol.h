#pragma once

#include "Expression.h"
#include "counter.h"
#include "utils/has_print_on.h"
#include <cstring>

namespace symbolic {
using utils::has_print_on::operator<<;

class Symbol : public ExpressionTag
{
 public:
  static constexpr precedence s_precedence = precedence::symbol;

 private:
  int const id_;
  char const* const name_;
  mutable int value_;

 public:
  template<auto Id = int{}, typename T = decltype([]{})>
  consteval Symbol(char const* name) : id_(metahack::unique_id<Id, T>()), name_(name), value_(0.0) { }

  consteval int id() const
  {
    return id_;
  }

  char const* name() const
  {
    return name_;
  }

  Symbol const& operator=(int value) const
  {
    value_ = value;
    return *this;
  }

#ifdef CWDEBUG
  bool needs_parens(before_or_after position, precedence prec) const { return prec < s_precedence; }
  bool needs_parens(precedence prec) const { return prec < s_precedence; }

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

template<>
struct is_symbol<Symbol> : std::true_type { };

template<typename E>
concept SymbolType = is_symbol_v<E>;

template<>
struct expression_order_less<Symbol, Symbol>
{
  Symbol const& symbol1_;
  Symbol const& symbol2_;

  consteval expression_order_less(Symbol const& symbol1, Symbol const& symbol2) : symbol1_(symbol1), symbol2_(symbol2) { }

  consteval bool operator()() const
  {
    return symbol1_.id() < symbol2_.id();
  }
};

} // namespace symbolic
