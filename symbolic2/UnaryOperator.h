#pragma once

#include "Expression.h"
#include "Hash.h"
#include "debug.h"

namespace symbolic2 {

class Sin;
class Cos;
class Atan;
class Log;

template<typename T>
concept UnaryOp = std::is_same_v<T, Sin> || std::is_same_v<T, Cos> || std::is_same_v<T, Atan> || std::is_same_v<T, Log>;

template<UnaryOp OpT>
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
};

} // namespace symbolic2

#include "Sin.h"
#include "Cos.h"
#include "Atan.h"
#include "Log.h"

namespace symbolic2 {

//static
template<UnaryOp OpT>
OpT const& UnaryOperator<OpT>::realize(Expression const& arg)
{
  DoutEntering(dc::symbolic, NAMESPACE_DEBUG::type_name_of<OpT>() << "::realize(\"" << arg << "\")");
  return static_cast<OpT const&>(get<OpT>(arg));
}

template<UnaryOp OpT>
bool UnaryOperator<OpT>::equals(Expression const& other) const
{
  OpT const* other_ = dynamic_cast<OpT const*>(&other);
  return other_ && arg_.equals(other_->arg_);
}

} // namespace symbolic2
