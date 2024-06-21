#pragma once

#include "utils/REMOVE_TRAILING_COMMA.h"
#include "debug.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

// Declare members of class Class##Style.
// These members are protected and begin with `m_`.
#define CAIROWINDOW_DECLARE_MEMBER(type, member, ...) \
  type m_##member;

// Declare members of struct Class##StyleParams.
// These members are public and have a default value.
// They can be changed using designated initialization.
#define CAIROWINDOW_DECLARE_STYLE_MEMBER(type, member, undval, Defaults, ...) type member = Defaults::member;

// Declare members of struct Class##StyleParamsDelta.
// These members are public and are initialized with a magic value meaning `undefined`.
// They can be changed using designated initialization.
#define CAIROWINDOW_DECLARE_STYLE_MEMBER_UNDEFINED(type, member, undefined_magic, ...) type member = undefined_magic;

// Used to build the initializer list of Class##Style classes.
#define CAIROWINDOW_INITIALIZER_LIST(type, member, ...) m_##member(params.member),

// Used to build the argument list of the base class of Class##Style classes that have a base class.
#define CAIROWINDOW_BASECLASS_PARAM_LIST(type, member, ...) params.member,

// Used to parse delta's in Class##Style classes.
#define CAIROWINDOW_UPDATE_STYLE_FROM_DELTA(type, member, undefined_magic, ...) \
  if (delta.member != undefined_magic) style.m_##member = delta.member;

// Declare member accessors of Class##Style classes.
#define CAIROWINDOW_DECLARE_MEMBER_ACCESSOR(type, member, ...) \
  type const& member() const { return m_##member; }

// Used to build the designated initialization list that copies from Defaults.
#define CAIROWINDOW_COPY_DEFAULTS(type, member, undval, Defaults, ...) \
  .member = Defaults::member,

// Used to build initializer list of Class##Style that copies from Defaults.
#define CAIROWINDOW_INITIALIZER_LIST_COPY_DEFAULTS(type, member, undval, Defaults, ...) \
  m_##member(Defaults::member),

#ifdef CWDEBUG
// Used for debug purposes; print members of the ostream `os`.
#define CAIROWINDOW_PRINT_ON(type, member, ...) os << #member ":" << member << "; ";
#define CAIROWINDOW_PRINT_m_ON(type, member, ...) os << #member ":" << m_##member << "; ";
#endif

// Macro used by DECLARE_STYLE and DECLARE_STYLE_WITH_BASE to declare the structs Class##StyleParams and Class##StyleParamsDelta.
#define CAIROWINDOW_DECLARE_STYLE_PARAMS(Class, Defaults) \
  using utils::has_print_on::operator<<; \
  struct Class##StyleParams { \
    cairowindow_##Class##_FOREACH_STYLE_MEMBER(CAIROWINDOW_DECLARE_STYLE_MEMBER, Defaults) \
    CWDEBUG_ONLY(void print_on(std::ostream& os) const { cairowindow_##Class##_FOREACH_STYLE_MEMBER(CAIROWINDOW_PRINT_ON) }) \
  }; \
  struct Class##StyleParamsDelta { \
    cairowindow_##Class##_FOREACH_STYLE_MEMBER(CAIROWINDOW_DECLARE_STYLE_MEMBER_UNDEFINED) \
  }

// Declare a new style class Class##Style.
//
// The macros cairowindow_##Class##_FOREACH_MEMBER(X, ...) and cairowindow_##Class##_FOREACH_STYLE_MEMBER(X, ...)
// must already have been defined before using this macro.
//
// Default must be a struct with static constexpr members of the same type and name,
// having the required default values.
#define DECLARE_STYLE(Class, Defaults) \
  CAIROWINDOW_DECLARE_STYLE_PARAMS(Class, Defaults); \
  using utils::has_print_on::operator<<; \
  class Class##Style { \
   protected: \
    cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_DECLARE_MEMBER) \
   public: \
    constexpr Class##Style(Class##StyleParams params) : \
      REMOVE_TRAILING_COMMA(cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_INITIALIZER_LIST)) { } \
    constexpr Class##Style() : Class##Style(Class##StyleParams{cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_COPY_DEFAULTS, Defaults)}) { } \
    Class##Style operator()(Class##StyleParamsDelta delta) const { \
      Class##Style style(*this); \
      cairowindow_##Class##_FOREACH_STYLE_MEMBER(CAIROWINDOW_UPDATE_STYLE_FROM_DELTA) \
      return style; \
    } \
    cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_DECLARE_MEMBER_ACCESSOR) \
    CWDEBUG_ONLY(void print_on(std::ostream& os) const { cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_PRINT_m_ON) }) \
  };

// Same, but derived from a base class.
// Base must already have been declared before using DECLARE_STYLE or DECLARE_STYLE_WITH_BASE.
#define DECLARE_STYLE_WITH_BASE(Class, Base, Defaults) \
  CAIROWINDOW_DECLARE_STYLE_PARAMS(Class, Defaults); \
  using utils::has_print_on::operator<<; \
  struct Class##StyleParamsPlus { \
    cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_DECLARE_STYLE_MEMBER, Defaults) \
  }; \
  class Class##Style : public Base##Style { \
   protected: \
    cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_DECLARE_MEMBER) \
   public: \
    constexpr Class##Style(Class##StyleParams params) : \
      Base##Style({REMOVE_TRAILING_COMMA(cairowindow_##Base##_FOREACH_STYLE_MEMBER(CAIROWINDOW_BASECLASS_PARAM_LIST))}) \
      REMOVE_TRAILING_COMMA(, cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_INITIALIZER_LIST)) { } \
    Class##Style() : Class##Style(Class##StyleParams{cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_COPY_DEFAULTS, Defaults)}) { } \
    Class##Style(Base##Style const& base) : Base##Style(base) \
      REMOVE_TRAILING_COMMA(, cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_INITIALIZER_LIST_COPY_DEFAULTS, Defaults)) { } \
    Class##Style(Base##Style base, Class##StyleParamsPlus params) : Base##Style(base) \
      REMOVE_TRAILING_COMMA(, cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_INITIALIZER_LIST)) { } \
    Class##Style operator()(Class##StyleParamsDelta delta) const \
    { \
      Class##Style style(*this); \
      cairowindow_##Class##_FOREACH_STYLE_MEMBER(CAIROWINDOW_UPDATE_STYLE_FROM_DELTA) \
      return style; \
    } \
    cairowindow_##Class##_FOREACH_MEMBER(CAIROWINDOW_DECLARE_MEMBER_ACCESSOR) \
    CWDEBUG_ONLY(void print_on(std::ostream& os) const { cairowindow_##Class##_FOREACH_STYLE_MEMBER(CAIROWINDOW_PRINT_m_ON) }) \
  };
