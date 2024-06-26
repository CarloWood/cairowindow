#pragma once

#include "Constant.h"

namespace symbolic {

class Symbol : public Expression
{
 private:
  std::string const name_;

  mutable double value_{};
  mutable uint64_t cached_hash_ = 0;

 private:
  template<typename T, typename... Args>
  friend Expression const& Expression::get(Args&&... args);
  Symbol(std::string const& name) : name_(name) { }

  Precedence precedence() const override final { return Precedence::symbol; }

 public:
  static Symbol const& realize(std::string const& name)
  {
    DoutEntering(dc::symbolic, "Symbol::realize(\"" << name << "\")");
    return static_cast<Symbol const&>(get<Symbol>(name));
  }

  ExpressionType type() const override final { return symbolT; }

  uint64_t hash() const override final
  {
    if (cached_hash_ == 0)
    {
      cached_hash_ = hash_v<Symbol>;
      boost::hash_combine(cached_hash_, std::hash<std::string>{}(name_));
    }

    return cached_hash_;
  }

  bool equals(Expression const& other) const override final;

  double evaluate() const override final
  {
    return value_;
  }

  Expression const& derivative(Symbol const& symbol) const override final
  {
    return &symbol == this ? Constant::s_cached_one : Constant::s_cached_zero;
  }

  std::string const& name() const { return name_; }

  Symbol const& operator=(double value) const
  {
    value_ = value;
    return *this;
  }

  bool operator<(Symbol const& other) const { return name_ < other.name_; }

#ifdef SYMBOLIC_PRINTING
 public:
  void print_on(std::ostream& os) const override final
  {
    os << name_;
  }
#endif
};

} // namespace symbolic
