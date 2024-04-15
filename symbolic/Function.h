#pragma once

#include "Symbol.h"
#include "Expression.h"
#include "Hash.h"
#include "Precedence.h"
#include <string>

namespace symbolic {

class Function : public Expression
{
 private:
  std::string const name_;
  Expression const& definition_;

  mutable uint64_t cached_hash_ = 0;

 private:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Function(std::string const& name, Expression const& definition) : name_(name), definition_(definition) { }

  Precedence precedence() const override final { return Precedence::symbol; }

 public:
  static Function const& realize(std::string const& name, Expression const& definition)
  {
    DoutEntering(dc::symbolic, "Function::realize(\"" << name << "\", " << definition <<  ")");
    return static_cast<Function const&>(get<Function>(name, definition));
  }

  ExpressionType type() const override final { return symbolT; }

  uint64_t hash() const override final
  {
    if (cached_hash_ == 0)
    {
      cached_hash_ = hash_v<Function>;
      boost::hash_combine(cached_hash_, std::hash<std::string>{}(name_));
      boost::hash_combine(cached_hash_, definition_.hash());
    }

    return cached_hash_;
  }

  bool equals(Expression const& other) const override final
  {
    Function const* other_function = dynamic_cast<Function const*>(&other);
    bool equal = other_function && name_ == other_function->name_;
    // Don't define two different function with the same name.
    ASSERT(!equal || definition_.equals(other_function->definition_));
    return equal;
  }

  double evaluate() const override final
  {
    return definition_.evaluate();
  }

  Function const& derivative(Symbol const& symbol) const override final
  {
    return realize("∂" + name_ + "/∂" + symbol.name(), definition_.derivative(symbol));
  }

  Expression const& definition() const { return definition_; }

#ifdef SYMBOLIC_PRINTING
 public:
  void print_on(std::ostream& os) const override final
  {
    os << name_;
  }
#endif
};

} // namespace symbolic
