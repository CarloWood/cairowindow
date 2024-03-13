#pragma once

#include "utils/macros.h"

namespace metahack {

template<auto Id>
struct counter
{
  using tag = counter;

  struct generator
  {
    friend consteval auto is_defined(tag) { return true; }
  };

PRAGMA_DIAGNOSTIC_PUSH_IGNORE_non_template_friend
  friend consteval auto is_defined(tag);
PRAGMA_DIAGNOSTIC_POP

  template<typename Tag = tag, auto = is_defined(Tag{})>
  static consteval bool exists(auto)
  {
    return true;
  }

  static consteval bool exists(...)
  {
    return generator(), false;
  }
};

template<auto Id = int{}, typename = decltype([]{})>
consteval auto unique_id()
{
  if constexpr (not counter<Id>::exists(Id))
    return Id;
  else
    return unique_id<Id + 1>();
}

} // namespace metahack
