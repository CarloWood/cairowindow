#pragma once

#include "Precedence.h"
#include <unordered_set>
#include <memory>
#include "debug.h"
#ifdef SYMBOLIC_PRINTING
#include "utils/has_print_on.h"
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
#endif

static constexpr int binary_op = 2;     //    v
static constexpr int unary_op = 4;      //  v

// A helper enum to determine the type and/or global order of available expression classes.
enum ExpressionType
{
  // Special cases:
  constantT =  0,                       // 00000
  symbolT = 1,                          // 00001
  // Binary operators:
  powerT = binary_op,                   // 00010
  productT = binary_op | 1,             // 00011
  // Unary operators:
  sinT = unary_op,                      // 00100
  cosT = unary_op | 1,                  // 00101
  logT = unary_op | 8,                  // 01100
  atanT = unary_op | 9,                 // 01101
  // Binary operators:
  sumT = binary_op | 16,                // 10010
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

  static void dump_database();

  static bool is_unary_op(ExpressionType type) { return (type & unary_op) != 0; }
  static bool is_binary_op(ExpressionType type) { return (type & binary_op) != 0; }

  bool is_unary_op() const { return (type() & unary_op) != 0; }
  bool is_binary_op() const { return (type() & binary_op) != 0; }

  bool is_constant() const { return type() == constantT; }
  bool is_symbol() const { return type() == symbolT; }
  bool is_power() const { return type() == powerT; }
  bool is_product() const { return type() == productT; }
  bool is_sum() const { return type() == sumT; }

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
