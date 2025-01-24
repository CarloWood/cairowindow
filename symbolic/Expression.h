#pragma once

#include "Precedence.h"
#include <unordered_set>
#include <memory>
#include "debug.h"
#ifdef SYMBOLIC_PRINTING
#include "utils/has_print_on.h"
#include "utils/iomanip.h"
#include <iostream>
#endif

#ifdef CWDEBUG
NAMESPACE_DEBUG_CHANNELS_START
extern channel_ct symbolic;
NAMESPACE_DEBUG_CHANNELS_END
#endif

namespace symbolic {
#ifdef SYMBOLIC_PRINTING
using utils::has_print_on::operator<<;

class UseUtf8 : public utils::iomanip::Unsticky<>
{
 private:
  static utils::iomanip::Index s_index;

 public:
  UseUtf8(long iword = 1L) : utils::iomanip::Unsticky<>(s_index, iword) { }

  static long get_iword_value(std::ostream& os) { return get_iword_from(os, s_index); }
};
#endif

static constexpr int binary_op = 8;     //   v
static constexpr int unary_op = 16;     //  v
static constexpr int binary_op_sum = binary_op | 32;
                                        // v v

// A helper enum to determine the type and/or global order of available expression classes.
enum ExpressionType
{                                       // .---- binary_op_sum
                                        // |.--- unary_op
                                        // ||.-- binary_op
                                        // vvv
  // Special cases:
  constantT = 0,                        // 000000
  symbolT   = 1,                        // 000001
  functionT = 2,                        // 000010
  // Binary operators:
  powerT   = binary_op,                 // 001000
  productT = binary_op | 1,             // 001001
  // Unary operators:
  expT  = unary_op,                     // 010000
  sinT  = unary_op | 1,                 // 010001
  cosT  = unary_op | 2,                 // 010010
  logT  = unary_op | 3,                 // 010011
  atanT = unary_op | 4,                 // 010100
  // Binary operators:
  sumT = binary_op_sum,                 // 101000
};

class Constant;
class Symbol;
struct KeyEqual;

// Everything is an Expression, even a Constant.
//
class Expression
{
 private:
  static std::unordered_set<std::unique_ptr<Expression>, std::hash<std::unique_ptr<Expression>>, KeyEqual> s_database_;

 public:
  virtual ~Expression() = default;

  virtual ExpressionType type() const = 0;
  virtual uint64_t hash() const = 0;
  virtual bool equals(Expression const& other) const = 0;
  virtual Expression const& arg1() const { ASSERT(false); AI_NEVER_REACHED }
  virtual Expression const& arg2() const { ASSERT(false); AI_NEVER_REACHED }
  virtual Expression const& get_nonconstant_factor() const { return *this; }
  virtual Constant const& get_constant_factor() const;
  virtual Expression const& get_base() const { return *this; }
  virtual Constant const& get_exponent() const;
  virtual Precedence precedence() const = 0;
  virtual double evaluate() const = 0;
  virtual Expression const& derivative(Symbol const& symbol) const = 0;
  virtual bool is_zero_function() const { return false; }
  virtual bool is_one_function() const { return false; }
  virtual bool is_minus_one_function() const { return false; }
  virtual Expression const* substitute(Expression const& replace, Expression const& with) const = 0;

  static void dump_database();

  static bool is_unary_op(ExpressionType type) { return (type & unary_op) != 0; }
  static bool is_binary_op(ExpressionType type) { return (type & binary_op) != 0; }

  bool is_unary_op() const { return (type() & unary_op) != 0; }
  bool is_binary_op() const { return (type() & binary_op) != 0; }

  bool is_constant() const { return type() == constantT || is_zero_function(); }
  bool is_symbol() const { return type() == symbolT; }
  bool is_function() const { return type() == functionT; }
  bool is_power() const { return type() == powerT; }
  bool is_product() const { return type() == productT; }
  bool is_sum() const { return type() == sumT; }

  Expression const& subs(Expression const& replace, Expression const& with) const
  {
    Expression const* new_expression = substitute(replace, with);
    return new_expression ? *new_expression : *this;
  }

 protected:
  template<typename T, typename... Args>
  static Expression const& get(Args&&... args);

#ifdef SYMBOLIC_PRINTING
 public:
  virtual void print_on(std::ostream& os) const = 0;
#endif
};

} // namespace symbolic

namespace std {

template<> struct hash<std::unique_ptr<symbolic::Expression>>
{
  uint64_t operator()(std::unique_ptr<symbolic::Expression> const& e) const
  {
    return e->hash();
  }
};

} // namespace std

namespace symbolic {

struct KeyEqual
{
  bool operator()(std::unique_ptr<Expression> const& key1, std::unique_ptr<Expression> const& key2) const
  {
    return key1->equals(*key2);
  }
};

template<typename T, typename... Args>
Expression const& Expression::get(Args&&... args)
{
  auto pib = s_database_.emplace(new T(std::forward<Args>(args)...));
//  if (pib.second)
//    Dout(dc::notice, "Inserted " << *pib.first->get() << " into the database. Hash = " << std::hash<std::unique_ptr<Expression>>{}(*pib.first));
  return *pib.first->get();
}

} // namespace symbolic
