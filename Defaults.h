#pragma once

#define DEFAULT_FROM(Defaults, member) \
  utils::defaults::default_##member##_from<Defaults, decltype(member)>()

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
