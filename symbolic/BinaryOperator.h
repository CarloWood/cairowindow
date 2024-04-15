#pragma once

#include "Symbol.h"
#include "Constant.h"
#include "Expression.h"
#include "Hash.h"

namespace symbolic {

class Sum;
class Product;
class Power;

template<typename T>
concept BinaryOp = std::is_same_v<T, Sum> || std::is_same_v<T, Product> || std::is_same_v<T, Power>;

template<BinaryOp OpT>
class BinaryOperator : public Expression
{
 protected:
  Expression const& arg1_;
  Expression const& arg2_;

  mutable uint64_t cached_hash_ = 0;

 protected:
  BinaryOperator(Expression const& arg1, Expression const& arg2) : arg1_(arg1), arg2_(arg2) { }

 public:
  static OpT const& realize(Expression const& arg1, Expression const& arg2);
  static bool is_less_exact(Expression const& arg1, Expression const& arg2);
  static bool is_less(Expression const& arg1, Expression const& arg2);

  uint64_t hash() const override final
  {
    if (cached_hash_ == 0)
    {
      cached_hash_ = hash_v<OpT>;
      boost::hash_combine(cached_hash_, arg1_.hash());
      boost::hash_combine(cached_hash_, arg2_.hash());
    }

    return cached_hash_;
  }

  Expression const& arg1() const override final { return arg1_; }
  Expression const& arg2() const override final { return arg2_; }

  bool equals(Expression const& other) const override final;
};

template<BinaryOp OpT>
bool BinaryOperator<OpT>::is_less_exact(Expression const& arg1, Expression const& arg2)
{
  if constexpr (std::is_same_v<OpT, Power>)
  {
    // This should never be called for an Power.
    ASSERT(false);
    return false;
  }
  else
  {
    ExpressionType type = arg1.type();
    ExpressionType type2 = arg2.type();

    if (type != type2)
    {
      if constexpr (std::is_same_v<OpT, Product>)
      {
        if (type2 == productT)
          return true;
        else if (type == productT)
          return false;
      }
      return type < type2;
    }

    if (type == constantT)
    {
      Constant const& constant1 = static_cast<Constant const&>(arg1);
      Constant const& constant2 = static_cast<Constant const&>(arg2);
      return constant1 < constant2;
    }
    else if (type == symbolT)
    {
      Symbol const& symbol1 = static_cast<Symbol const&>(arg1);
      Symbol const& symbol2 = static_cast<Symbol const&>(arg2);
      return symbol1 < symbol2;
    }
    else if (is_binary_op(type))
    {
      return is_less_exact(arg1.arg1(), arg2.arg1()) || (!is_less_exact(arg2.arg1(), arg1.arg1()) && is_less_exact(arg1.arg2(), arg2.arg2()));
    }
    else // unary operator
    {
      return is_less_exact(arg1.arg1(), arg2.arg1());
    }
  }
}

template<BinaryOp OpT>
bool BinaryOperator<OpT>::is_less(Expression const& arg1, Expression const& arg2)
{
  static_assert(!std::is_same_v<OpT, Power>, "we should never get here?!");

  if constexpr (std::is_same_v<OpT, Sum>)
    return is_less_exact(arg1.get_nonconstant_factor(), arg2.get_nonconstant_factor());
  else if (arg1.is_constant() && arg2.is_constant())
    return false;
  else
    return is_less_exact(arg1.get_base(), arg2.get_base());
}

} // namespace symbolic

#include "Power.h"

namespace symbolic {

//static
template<BinaryOp OpT>
OpT const& BinaryOperator<OpT>::realize(Expression const& arg1, Expression const& arg2)
{
  DoutEntering(dc::symbolic|continued_cf, NAMESPACE_DEBUG::type_name_of<OpT>() << "::realize(" << arg1 << ", " << arg2 << ") --> ");
  OpT const& result = static_cast<OpT const&>(get<OpT>(arg1, arg2));
  Dout(dc::finish, result);
  return result;
}

template<BinaryOp OpT>
bool BinaryOperator<OpT>::equals(Expression const& other) const
{
  OpT const* other_ = dynamic_cast<OpT const*>(&other);
  return other_ && arg1_.equals(other_->arg1_) && arg2_.equals(other_->arg2_);
}

} // namespace symbolic
