#pragma once

#include "Constant.h"

namespace symbolic2 {

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
  Symbol(std::string const& name, double value) : name_(name), value_(value) { }

 public:
  static Symbol const& realize(std::string const& name) { return static_cast<Symbol const&>(get<Symbol>(name)); }
  static Symbol const& realize(std::string const& name, double value) { return static_cast<Symbol const&>(get<Symbol>(name, value)); }

  uint64_t hash() const override
  {
    if (cached_hash_ == 0)
    {
      cached_hash_ = hash_v<Symbol>;
      boost::hash_combine(cached_hash_, std::hash<std::string>{}(name_));
    }

    return cached_hash_;
  }

  bool equals(Expression const& other) const override;

  double evaluate() const override
  {
    return value_;
  }

  Expression const& differentiate(Symbol const& symbol) const override
  {
    return &symbol == this ? Constant::s_cached_one : Constant::s_cached_zero;
  }

  Symbol const& operator=(double value) const
  {
    value_ = value;
    return *this;
  }

#ifdef SYMBOLIC2_PRINTING
 public:
  void print_on(std::ostream& os) const override
  {
    os << name_;
  }
#endif
};

} // namespace symbolic2
