#pragma once

#include "Expression.h"
#include "counter.h"
#include "IdRange.h"
#include "precedence.h"
#ifdef SYMBOLIC_PRINTING
#include "utils/has_print_on.h"
#include <iostream>
#endif

namespace symbolic {
#ifdef SYMBOLIC_PRINTING
using utils::has_print_on::operator<<;
#endif

template<int Id, typename T>
inline consteval auto make_symbol();

template<int Id, typename T>
inline auto make_symbol(char const* name);

class SymbolRegistry
{
 protected:
  static void register_symbol(int id, char const* name);
  static char const* get_name(int id);
  static void set_value(int id, double value);

 public:
  static double get_value(int id);
};

template<int Id>
class Symbol : public ExpressionTag, public SymbolRegistry
{
 public:
  static constexpr precedence s_precedence = precedence::symbol;
  static constexpr int s_id = Id;
  static constexpr IdRange<Id, Id + 1> id_range{};

 private:
  // The constructor is private. Use make_symbol to create a new symbol.

  template<int Id2, typename T2>
  friend consteval auto make_symbol();

  template<int Id2, typename T2>
  friend auto make_symbol(char const* name);

  Symbol(char const* name) { SymbolRegistry::register_symbol(s_id, name); }

  // Or, use these:
  consteval Symbol() = default;
 public:
  static void register_name(char const* name) { SymbolRegistry::register_symbol(s_id, name); }

 public:
  static char const* name() { return SymbolRegistry::get_name(s_id); }

  Symbol const& operator=(double value) const
  {
    SymbolRegistry::set_value(s_id, value);
    return *this;
  }

  static constexpr Symbol instance() { return {}; }

#ifdef SYMBOLIC_PRINTING
  static void print_on(std::ostream& os)
  {
    os << name();
  }
#endif
};

template<int Id = 0, typename T = decltype([]{})>
auto make_symbol(char const* name)
{
  auto symbol = Symbol<metahack::unique_id<Id, T>()>(name);
  return symbol;
}

template<int Id = 0, typename T = decltype([]{})>
consteval auto make_symbol()
{
  auto symbol = Symbol<metahack::unique_id<Id, T>()>();
  return symbol;
}

template<int Id1, int Id2>
constexpr bool is_same_expression(Symbol<Id1> const& arg1, Symbol<Id2> const& arg2) { return arg1.id() == arg2.id(); }

} // namespace symbolic
