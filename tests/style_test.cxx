#include "sys.h"
#include "cairowindow/Style.h"
#include "utils/has_print_on.h"
#include "utils/REMOVE_TRAILING_COMMA.h"
#include "debug.h"

static constexpr int magic_undefined1 = 0x123456;
static constexpr int magic_undefined2 = -1;
static constexpr int magic_undefined3 = -2;

// List the members of BarStyle.
#define cairowindow_Bar_FOREACH_MEMBER(X, ...) \
  X(int, A, magic_undefined2, __VA_ARGS__) \
  X(int, B, magic_undefined3, __VA_ARGS__)

// Mandatory macro.
#define cairowindow_Bar_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Bar_FOREACH_MEMBER(X, __VA_ARGS__)

// List the additional members of FooStyle.
#define cairowindow_Foo_FOREACH_MEMBER(X, ...) \
  X(int, C, magic_undefined1, __VA_ARGS__) \
  X(int, D, magic_undefined3, __VA_ARGS__)

// FooStyle is derived from BarStyle.
#define cairowindow_Foo_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Bar_FOREACH_MEMBER(X, __VA_ARGS__) \
  cairowindow_Foo_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for BarStyle.
struct BarStyleParamsDefault
{
  static int const A = 101;
  static int const B = 102;
};

// Define default values for FooStyle.
struct FooStyleParamsDefault : BarStyleParamsDefault
{
  static int const B = 202;   // Change default from BarStyleParamsDefault.
  static int const C = 203;
  static int const D = 204;
};

// Declare BarStyle.
DECLARE_STYLE(Bar, BarStyleParamsDefault);

// Declare FooStyle, derived from BarStyle.
DECLARE_STYLE_WITH_BASE(Foo, Bar, FooStyleParamsDefault);

// A test function taking a FooStyle.
void f(FooStyle const& style)
{
  Dout(dc::notice, "f(" << style << ")");
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());
  Dout(dc::notice, "Entering main()");

  // Construct a FooStyle, using the default values from FooStyleParamsDefault,
  // except for A and C that are overwritten using designated initialization.
  FooStyle foo_style({.A = 3, .C = 42});

  // Pass foo_style to f(), but make another change using designated initialization.
  f(foo_style({.B = 13, .C = 2}));

  Dout(dc::notice, "foo_style = " << foo_style);

  Dout(dc::notice, "Leaving main()");
}
