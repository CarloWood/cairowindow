#pragma once

#include "Symbol.h"
#include "Expression.h"
#include "Hash.h"
#include "Precedence.h"
#include <string>
#include <optional>
#ifdef SYMBOLIC_PRINTING
#include "iomanip_fulldef.h"
#include <sstream>
#endif

namespace symbolic {

class Function : public Expression
{
 public:
  static constexpr ExpressionType expression_type = functionT;

 private:
  std::string const name_;
  Expression const& definition_;

  mutable uint64_t cached_hash_ = 0;
  mutable std::optional<double> evaluation_;
  mutable Function const* derivative_cache_ = nullptr;

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

  ExpressionType type() const override final { return expression_type; }

  uint64_t hash() const override final
  {
    if (cached_hash_ == 0)
    {
      cached_hash_ = hash_v<Function>;
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

  void reset_evaluation() const
  {
    evaluation_.reset();
    if (derivative_cache_)
      derivative_cache_->reset_evaluation();
  }

  double evaluate() const override final
  {
    if (!evaluation_.has_value())
      evaluation_ = definition_.evaluate();
    return *evaluation_;
  }

  Function const& derivative(Symbol const& symbol) const override final
  {
    if (!derivative_cache_)
      derivative_cache_ = &realize("∂" + name_ + "/∂" + symbol.name(), definition_.derivative(symbol));

    return *derivative_cache_;
  }

  Expression const& definition() const { return definition_; }

  bool operator<(Function const& other) const { return hash() < other.hash(); }

#ifdef SYMBOLIC_PRINTING
 public:
  std::string to_string() const
  {
    std::ostringstream oss;
    print_on(oss << fulldef);
    return oss.str();
  }

  void print_on(std::ostream& os) const override final
  {
    if (IOManipFullDef::is_full_def(os))
      os << '(' << definition_ << ')';
    else
      os << name_;
  }
#endif
};

} // namespace symbolic
