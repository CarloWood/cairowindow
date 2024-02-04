#pragma once

#define DECLARE_WITH_DEFAULT_FROM(type, member, Defaults) \
  type member = utils::defaults::default_##member##_from<Defaults, type>()

#define DECLARE_DEFAULTS_HAS_MEMBER(member) \
  namespace utils::defaults { \
    template<typename Defaults, typename T> \
    static constexpr T default_##member##_from() \
    { \
      if constexpr (requires { Defaults::member; }) \
        return Defaults::member; \
      else \
        return {}; \
    } \
  } // namespace utils::defaults
