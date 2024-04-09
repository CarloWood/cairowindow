#pragma once

#include <unordered_set>
#include <memory>
#include "debug.h"
#ifdef SYMBOLIC2_PRINTING
#include "utils/has_print_on.h"
#include <iostream>
#endif

namespace symbolic2 {
#ifdef SYMBOLIC2_PRINTING
using utils::has_print_on::operator<<;
#endif

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

  virtual uint64_t hash() const = 0;
  virtual bool equals(Expression const& other) const = 0;
  virtual double evaluate() const = 0;
  virtual Expression const& differentiate(Symbol const& symbol) const = 0;

  static void dump_database();

 protected:
  template<typename T, typename... Args>
  static Expression const& get(Args&&... args);

#ifdef SYMBOLIC2_PRINTING
 public:
  virtual void print_on(std::ostream& os) const  = 0;
#endif
};

} // namespace symbolic2

namespace std {

template<> struct hash<std::unique_ptr<symbolic2::Expression>>
{
  uint64_t operator()(std::unique_ptr<symbolic2::Expression> const& e) const
  {
    return e->hash();
  }
};

} // namespace std

namespace symbolic2 {

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

} // namespace symbolic2
