#pragma once

#include "Expression.h"
#include "Hash.h"
#include <boost/functional/hash.hpp>
#include "debug.h"

namespace symbolic {

template<typename OpT>
class UnaryOperator : public Expression
{
 protected:
  Expression const& arg_;

  mutable uint64_t cached_hash_ = 0;

 protected:
  UnaryOperator(Expression const& arg) : arg_(arg) { }

 public:
  static OpT const& realize(Expression const& arg);

  uint64_t hash() const override final
  {
    if (cached_hash_ == 0)
    {
      cached_hash_ = hash_v<OpT>;
      boost::hash_combine(cached_hash_, arg_.hash());
    }

    return cached_hash_;
  }

  Expression const& arg1() const override final { return arg_; }

  bool equals(Expression const& other) const override final;

  Expression const* substitute(Expression const& replace, Expression const& with) const override final
  {
    Expression const* new_arg = arg_.substitute(replace, with);
    return new_arg ? &realize(*new_arg) : nullptr;
  }
};

} // namespace symbolic

namespace symbolic {

//static
template<typename OpT>
OpT const& UnaryOperator<OpT>::realize(Expression const& arg)
{
  DoutEntering(dc::symbolic, NAMESPACE_DEBUG::type_name_of<OpT>() << "::realize(\"" << arg << "\")");
  return static_cast<OpT const&>(get<OpT>(arg));
}

template<typename OpT>
bool UnaryOperator<OpT>::equals(Expression const& other) const
{
  OpT const* other_ = dynamic_cast<OpT const*>(&other);
  return other_ && arg_.equals(other_->arg_);
}

} // namespace symbolic
