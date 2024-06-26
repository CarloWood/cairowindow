#include "sys.h"
#include "symbolic.h"
#include "utils/square.h"
#include "debug.h"

using namespace symbolic;

void test_print(char const* what, Expression const& expression, std::string expect_string, int line)
{
#ifdef SYMBOLIC_PRINTING
  std::ostringstream result;
  result << expression;
#ifdef CWDEBUG
  Dout(dc::notice, what << " = " << result.str());
  if (result.str() != expect_string)
  {
    DoutFatal(dc::core, "symbolic_testsuite.cxx:" << line << ": expected \"" << expect_string << "\" got \"" << result.str() << "\".");
  }
#else
  std::cout << what << " = " << result.str() << std::endl;
  assert(result.str() == expect_string);
#endif
#endif  // SYMBOLIC_PRINTING
}

void test_same_expression(Expression const& expression, Expression const& expect, int line)
{
  if (&expression != &expect)
  {
#ifdef SYMBOLIC_PRINTING
    DoutFatal(dc::core, "symbolic_testsuite.cxx:" << line << ": expected \"" << expression << "\" and \"" << expect << "\" to be the same expression.");
#else
    DoutFatal(dc::core, "symbolic_testsuite.cxx:" << line << ": unexpected result string.");
#endif
  }
}

#define TEST(what, expression, expect_string) \
  test_print((what), (expression), (expect_string), __LINE__)

#define TESTT(expression, expect) \
  test_same_expression((expression), (expect), __LINE__)

#define TESTS(x, y) TEST(#x, x, y)

void compare_equal1(Expression const& expression)
{
  // "Expected E1 == E1";
  ASSERT(!Sum::is_less(expression, expression));
}

enum What
{
  ConstantsCompareNotEqual = 0,
  ConstantsCompareEqual,
  ConstantsCompareEitherWay
};

void test_compare_less(Expression const& E1, Expression const& E2, int line, What constants_compare_equal = ConstantsCompareEitherWay)
{
  bool fail = false;

  if (constants_compare_equal == ConstantsCompareNotEqual || constants_compare_equal == ConstantsCompareEitherWay)
  {
    bool should_be_true = Sum::is_less_exact(E1, E2);
    if (!should_be_true)
    {
      Dout(dc::warning, "symbolic_testsuite.cxx:" << line << ": expected \"" << E1 << "\" to be less than \"" << E2 << "\" when constants do not compare equal");
      fail = true;
    }
    bool should_be_false = Sum::is_less_exact(E2, E1);
    if (should_be_false)
    {
      Dout(dc::warning, "symbolic_testsuite.cxx:" << line << ": expected \"" << E2 << "\" NOT to be less than \"" << E1 << "\" when constants do not compare equal");
      fail = true;
    }
  }
  if (constants_compare_equal == ConstantsCompareEqual || constants_compare_equal == ConstantsCompareEitherWay)
  {
    bool should_be_true = Sum::is_less(E1, E2);
    if (!should_be_true)
    {
      Dout(dc::warning, "symbolic_testsuite.cxx:" << line << ": expected \"" << E1 << "\" to be less than \"" << E2 << "\" when constants compare equal");
      fail = true;
    }
    bool should_be_false = Sum::is_less(E2, E1);
    if (should_be_false)
    {
      Dout(dc::warning, "symbolic_testsuite.cxx:" << line << ": expected \"" << E2 << "\" NOT to be less than \"" << E1 << "\" when constants compare equal");
      fail = true;
    }
  }
  ASSERT(!fail);
}

#define compare_less(E1, E2) test_compare_less(E1, E2, __LINE__)
#define compare_less_ConstantsCompareEqual(E1, E2) test_compare_less(E1, E2, __LINE__, ConstantsCompareEqual)
#define compare_less_ConstantsCompareNotEqual(E1, E2) test_compare_less(E1, E2, __LINE__, ConstantsCompareNotEqual)

void test_compare_equal(Expression const& E1, Expression const& E2, int line, What constants_compare_equal = ConstantsCompareEitherWay)
{
  bool fail = false;

  if (constants_compare_equal == ConstantsCompareNotEqual || constants_compare_equal == ConstantsCompareEitherWay)
  {
    bool should_be_false = Sum::is_less_exact(E1, E2);
    if (should_be_false)
    {
      Dout(dc::warning, "symbolic_testsuite.cxx:" << line << ": expected \"" << E1 << "\" NOT to be less than \"" << E2 << "\" when constants do not compare equal");
      fail = true;
    }
    should_be_false = Sum::is_less_exact(E2, E1);
    if (should_be_false)
    {
      Dout(dc::warning, "symbolic_testsuite.cxx:" << line << ": expected \"" << E2 << "\" NOT to be less than \"" << E1 << "\" when constants do not compare equal");
      fail = true;
    }
  }
  if (constants_compare_equal == ConstantsCompareEqual || constants_compare_equal == ConstantsCompareEitherWay)
  {
    bool should_be_false = Sum::is_less(E1, E2);
    if (should_be_false)
    {
      Dout(dc::warning, "symbolic_testsuite.cxx:" << line << ": expected \"" << E1 << "\" NOT to be less than \"" << E2 << "\" when constants compare equal");
      fail = true;
    }
    should_be_false = Sum::is_less(E2, E1);
    if (should_be_false)
    {
      Dout(dc::warning, "symbolic_testsuite.cxx:" << line << ": expected \"" << E2 << "\" NOT to be less than \"" << E1 << "\" when constants compare equal");
      fail = true;
    }
  }
  ASSERT(!fail);
}

#define compare_equal(E1, E2) test_compare_equal(E1, E2, __LINE__)
#define compare_equal_ConstantsCompareEqual(E1, E2) test_compare_equal(E1, E2, __LINE__, ConstantsCompareEqual)
#define compare_equal_ConstantsCompareNotEqual(E1, E2) test_compare_equal(E1, E2, __LINE__, ConstantsCompareNotEqual)

void test_derivative(Expression const& expression, Symbol const& symbol)
{
  Dout(dc::notice, "∂/∂" << symbol << " (" << expression << ") = " << expression.derivative(symbol));
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  Constant const& zero = Constant::realize(0);
  Constant const& one = Constant::realize(1);
  Constant const& two = Constant::realize(2);
  Constant const& three = Constant::realize(3);
  Constant const& three_halfs = Constant::realize(3, 2);

  Symbol const& a = Symbol::realize("a");
  Symbol const& b = Symbol::realize("b");
  Symbol const& c = Symbol::realize("c");
  Symbol const& d = Symbol::realize("d");
  Symbol const& e = Symbol::realize("e");
  Symbol const& f = Symbol::realize("f");
  Symbol const& g = Symbol::realize("g");
  Symbol const& h = Symbol::realize("h");

  Symbol const& v = Symbol::realize("v");
  Symbol const& w = Symbol::realize("w");
  Symbol const& x = Symbol::realize("x");
  Symbol const& y = Symbol::realize("y");
  Symbol const& z = Symbol::realize("z");

  // Constant creation.
  auto constant = [](int e, int d = 1){ return Constant::realize(e, d); };

  TEST("42 / 1", constant(42, 1), "42");
  TEST("-42 / 1", constant(-42, 1), "-42");
  TEST("42 / -1", constant(42, -1), "-42");
  TEST("369 / 3", constant(369, 3), "123");
  TEST("-369 / -3", constant(-369, -3), "123");
  TEST("(2 * 3 * 7 * 13 / (2 * 7 * 11 * 17)", constant(2 * 3 * 7 * 13, 2 * 7 * 11 * 17), "39/187");

  // Constant negation.
  TEST("-(13/42)", -constant(13, 42), "-13/42");
  TEST("-(-13/42)", -constant(-13, 42), "13/42");

  Constant const& A = constant(2, 3 * 7);
  Constant const& B = constant(-5, 7 * 11);
  Constant const& C = constant(5, 7 * 11);

  // Addition of two Constant's.
  TEST("A(2/21) + B(-5/77)", A + B, "1/33");

  // Subtraction of two Constant's.
  TEST("A(2/21) - B(-5/77)", A - B, "37/231");
  TEST("A(2/21) - C(5/77)", A - C, "1/33");

  // Multiplication of two Constant's.
  TEST("A(2/21) * B(-5/77)", A * B, "-10/1617");

  // Division of two Constant's.
  TEST("A(2/21) / B(-5/77)", A / B, "-22/15");

  Constant const& C2 = constant(2);
  Constant const& C29 = constant(29);

  // Exponentiation of a Constant.
  TEST("B(-5/77) ^ 3", B^3, "-125/456533");
  TEST("C2(2) ^ C29(29)", C2^C29, "536870912");

  // Negation of a Constant.
  TEST("-B(-5/77)", -B, "5/77");

  // Negation of a Symbol.
  TEST("-x", -x, "-x");
  // Negation of a Negation.
  TEST("-(-x)", -(-x), "x");

  // Negation tests involving products.

  auto& minus_x = -x;
  TESTT(-x, minus_x);
  TESTT(-minus_x, x);

  auto& x_times_y = x * y;
  auto& minus__x_times_y = minus_x * y;
  auto& minus_y = -y;
  TESTT(x * minus_y, minus__x_times_y);

  // Multiplication of a symbol and a constant.
  TESTS(x * zero, "0");
  TESTS(zero * x, "0");
  TESTS(x * one, "x");
  TESTS(one * x, "x");
  TESTS((x^three_halfs) * zero, "0");
  TESTS(zero * (x^three_halfs), "0");
  TESTS((x^three_halfs) * one, "x^(3/2)");
  TESTS(one * (x^three_halfs), "x^(3/2)");

  // Multiplication of two symbols.
  TESTS(x * y, "x * y");
  TESTS(y * x, "x * y");
  TESTS(x * x, "x^2");

  // Multiplication of three symbols.
  char const* expected3 = "x * y * z";
  TESTS((x * y) * z, expected3);
  TESTS((x * z) * y, expected3);
  TESTS((y * x) * z, expected3);
  TESTS((y * z) * x, expected3);
  TESTS((z * x) * y, expected3);
  TESTS((z * y) * x, expected3);
  TESTS(x * (y * z), expected3);
  TESTS(x * (z * y), expected3);
  TESTS(y * (x * z), expected3);
  TESTS(y * (z * x), expected3);
  TESTS(z * (x * y), expected3);
  TESTS(z * (y * x), expected3);

  // Multiplication of four symbols.
  char const* expected4 = "w * x * y * z";
  TESTS(((w * x) * y) * z, expected4);
  TESTS((w * x) * (y * z), expected4);
  TESTS((w * (x * y)) * z, expected4);
  TESTS(w * ((x * y) * z), expected4);
  TESTS(w * (x * (y * z)), expected4);
  TESTS(((w * x) * z) * y, expected4);
  TESTS((w * x) * (z * y), expected4);
  TESTS((w * (x * z)) * y, expected4);
  TESTS(w * ((x * z) * y), expected4);
  TESTS(w * (x * (z * y)), expected4);
  TESTS(((w * y) * x) * z, expected4);
  TESTS((w * y) * (x * z), expected4);
  TESTS((w * (y * x)) * z, expected4);
  TESTS(w * ((y * x) * z), expected4);
  TESTS(w * (y * (x * z)), expected4);
  TESTS(((w * y) * z) * x, expected4);
  TESTS((w * y) * (z * x), expected4);
  TESTS((w * (y * z)) * x, expected4);
  TESTS(w * ((y * z) * x), expected4);
  TESTS(w * (y * (z * x)), expected4);
  TESTS(((w * z) * x) * y, expected4);
  TESTS((w * z) * (x * y), expected4);
  TESTS((w * (z * x)) * y, expected4);
  TESTS(w * ((z * x) * y), expected4);
  TESTS(w * (z * (x * y)), expected4);
  TESTS(((w * z) * y) * x, expected4);
  TESTS((w * z) * (y * x), expected4);
  TESTS((w * (z * y)) * x, expected4);
  TESTS(w * ((z * y) * x), expected4);
  TESTS(w * (z * (y * x)), expected4);
  TESTS(((x * w) * y) * z, expected4);
  TESTS((x * w) * (y * z), expected4);
  TESTS((x * (w * y)) * z, expected4);
  TESTS(x * ((w * y) * z), expected4);
  TESTS(x * (w * (y * z)), expected4);
  TESTS(((x * w) * z) * y, expected4);
  TESTS((x * w) * (z * y), expected4);
  TESTS((x * (w * z)) * y, expected4);
  TESTS(x * ((w * z) * y), expected4);
  TESTS(x * (w * (z * y)), expected4);
  TESTS(((x * y) * w) * z, expected4);
  TESTS((x * y) * (w * z), expected4);
  TESTS((x * (y * w)) * z, expected4);
  TESTS(x * ((y * w) * z), expected4);
  TESTS(x * (y * (w * z)), expected4);
  TESTS(((x * y) * z) * w, expected4);
  TESTS((x * y) * (z * w), expected4);
  TESTS((x * (y * z)) * w, expected4);
  TESTS(x * ((y * z) * w), expected4);
  TESTS(x * (y * (z * w)), expected4);
  TESTS(((x * z) * w) * y, expected4);
  TESTS((x * z) * (w * y), expected4);
  TESTS((x * (z * w)) * y, expected4);
  TESTS(x * ((z * w) * y), expected4);
  TESTS(x * (z * (w * y)), expected4);
  TESTS(((x * z) * y) * w, expected4);
  TESTS((x * z) * (y * w), expected4);
  TESTS((x * (z * y)) * w, expected4);
  TESTS(x * ((z * y) * w), expected4);
  TESTS(x * (z * (y * w)), expected4);
  TESTS(((y * w) * x) * z, expected4);
  TESTS((y * w) * (x * z), expected4);
  TESTS((y * (w * x)) * z, expected4);
  TESTS(y * ((w * x) * z), expected4);
  TESTS(y * (w * (x * z)), expected4);
  TESTS(((y * w) * z) * x, expected4);
  TESTS((y * w) * (z * x), expected4);
  TESTS((y * (w * z)) * x, expected4);
  TESTS(y * ((w * z) * x), expected4);
  TESTS(y * (w * (z * x)), expected4);
  TESTS(((y * x) * w) * z, expected4);
  TESTS((y * x) * (w * z), expected4);
  TESTS((y * (x * w)) * z, expected4);
  TESTS(y * ((x * w) * z), expected4);
  TESTS(y * (x * (w * z)), expected4);
  TESTS(((y * x) * z) * w, expected4);
  TESTS((y * x) * (z * w), expected4);
  TESTS((y * (x * z)) * w, expected4);
  TESTS(y * ((x * z) * w), expected4);
  TESTS(y * (x * (z * w)), expected4);
  TESTS(((y * z) * w) * x, expected4);
  TESTS((y * z) * (w * x), expected4);
  TESTS((y * (z * w)) * x, expected4);
  TESTS(y * ((z * w) * x), expected4);
  TESTS(y * (z * (w * x)), expected4);
  TESTS(((y * z) * x) * w, expected4);
  TESTS((y * z) * (x * w), expected4);
  TESTS((y * (z * x)) * w, expected4);
  TESTS(y * ((z * x) * w), expected4);
  TESTS(y * (z * (x * w)), expected4);
  TESTS(((z * w) * x) * y, expected4);
  TESTS((z * w) * (x * y), expected4);
  TESTS((z * (w * x)) * y, expected4);
  TESTS(z * ((w * x) * y), expected4);
  TESTS(z * (w * (x * y)), expected4);
  TESTS(((z * w) * y) * x, expected4);
  TESTS((z * w) * (y * x), expected4);
  TESTS((z * (w * y)) * x, expected4);
  TESTS(z * ((w * y) * x), expected4);
  TESTS(z * (w * (y * x)), expected4);
  TESTS(((z * x) * w) * y, expected4);
  TESTS((z * x) * (w * y), expected4);
  TESTS((z * (x * w)) * y, expected4);
  TESTS(z * ((x * w) * y), expected4);
  TESTS(z * (x * (w * y)), expected4);
  TESTS(((z * x) * y) * w, expected4);
  TESTS((z * x) * (y * w), expected4);
  TESTS((z * (x * y)) * w, expected4);
  TESTS(z * ((x * y) * w), expected4);
  TESTS(z * (x * (y * w)), expected4);
  TESTS(((z * y) * w) * x, expected4);
  TESTS((z * y) * (w * x), expected4);
  TESTS((z * (y * w)) * x, expected4);
  TESTS(z * ((y * w) * x), expected4);
  TESTS(z * (y * (w * x)), expected4);
  TESTS(((z * y) * x) * w, expected4);
  TESTS((z * y) * (x * w), expected4);
  TESTS((z * (y * x)) * w, expected4);
  TESTS(z * ((y * x) * w), expected4);
  TESTS(z * (y * (x * w)), expected4);

  // Test sorting of 8 symbols during multiplication.
  Expression const& expected8 = a * (b * (c * (d * (e * (f * (g * h))))));
  TESTT(a * (b * c * d * e * f * g * h), expected8);
  TESTT(b * (a * c * d * e * f * g * h), expected8);
  TESTT(c * (a * b * d * e * f * g * h), expected8);
  TESTT(d * (a * b * c * e * f * g * h), expected8);
  TESTT(e * (a * b * c * d * f * g * h), expected8);
  TESTT(f * (a * b * c * d * e * g * h), expected8);
  TESTT(g * (a * b * c * d * e * f * h), expected8);
  TESTT(h * (a * b * c * d * e * f * g), expected8);
  TESTT((a * b) * (c * d * e * f * g * h), expected8);
  TESTT((a * c) * (b * d * e * f * g * h), expected8);
  TESTT((a * d) * (b * c * e * f * g * h), expected8);
  TESTT((a * e) * (b * c * d * f * g * h), expected8);
  TESTT((a * f) * (b * c * d * e * g * h), expected8);
  TESTT((a * g) * (b * c * d * e * f * h), expected8);
  TESTT((a * h) * (b * c * d * e * f * g), expected8);
  TESTT((b * c) * (a * d * e * f * g * h), expected8);
  TESTT((b * d) * (a * c * e * f * g * h), expected8);
  TESTT((b * e) * (a * c * d * f * g * h), expected8);
  TESTT((b * f) * (a * c * d * e * g * h), expected8);
  TESTT((b * g) * (a * c * d * e * f * h), expected8);
  TESTT((b * h) * (a * c * d * e * f * g), expected8);
  TESTT((c * d) * (a * b * e * f * g * h), expected8);
  TESTT((c * e) * (a * b * d * f * g * h), expected8);
  TESTT((c * f) * (a * b * d * e * g * h), expected8);
  TESTT((c * g) * (a * b * d * e * f * h), expected8);
  TESTT((c * h) * (a * b * d * e * f * g), expected8);
  TESTT((d * e) * (a * b * c * f * g * h), expected8);
  TESTT((d * f) * (a * b * c * e * g * h), expected8);
  TESTT((d * g) * (a * b * c * e * f * h), expected8);
  TESTT((d * h) * (a * b * c * e * f * g), expected8);
  TESTT((e * f) * (a * b * c * d * g * h), expected8);
  TESTT((e * g) * (a * b * c * d * f * h), expected8);
  TESTT((e * h) * (a * b * c * d * f * g), expected8);
  TESTT((f * g) * (a * b * c * d * e * h), expected8);
  TESTT((f * h) * (a * b * c * d * e * g), expected8);
  TESTT((g * h) * (a * b * c * d * e * f), expected8);
  TESTT((a * b * c) * (d * e * f * g * h), expected8);
  TESTT((a * b * d) * (c * e * f * g * h), expected8);
  TESTT((a * b * e) * (c * d * f * g * h), expected8);
  TESTT((a * b * f) * (c * d * e * g * h), expected8);
  TESTT((a * b * g) * (c * d * e * f * h), expected8);
  TESTT((a * b * h) * (c * d * e * f * g), expected8);
  TESTT((a * c * d) * (b * e * f * g * h), expected8);
  TESTT((a * c * e) * (b * d * f * g * h), expected8);
  TESTT((a * c * f) * (b * d * e * g * h), expected8);
  TESTT((a * c * g) * (b * d * e * f * h), expected8);
  TESTT((a * c * h) * (b * d * e * f * g), expected8);
  TESTT((a * d * e) * (b * c * f * g * h), expected8);
  TESTT((a * d * f) * (b * c * e * g * h), expected8);
  TESTT((a * d * g) * (b * c * e * f * h), expected8);
  TESTT((a * d * h) * (b * c * e * f * g), expected8);
  TESTT((a * e * f) * (b * c * d * g * h), expected8);
  TESTT((a * e * g) * (b * c * d * f * h), expected8);
  TESTT((a * e * h) * (b * c * d * f * g), expected8);
  TESTT((a * f * g) * (b * c * d * e * h), expected8);
  TESTT((a * f * h) * (b * c * d * e * g), expected8);
  TESTT((a * g * h) * (b * c * d * e * f), expected8);
  TESTT((b * c * d) * (a * e * f * g * h), expected8);
  TESTT((b * c * e) * (a * d * f * g * h), expected8);
  TESTT((b * c * f) * (a * d * e * g * h), expected8);
  TESTT((b * c * g) * (a * d * e * f * h), expected8);
  TESTT((b * c * h) * (a * d * e * f * g), expected8);
  TESTT((b * d * e) * (a * c * f * g * h), expected8);
  TESTT((b * d * f) * (a * c * e * g * h), expected8);
  TESTT((b * d * g) * (a * c * e * f * h), expected8);
  TESTT((b * d * h) * (a * c * e * f * g), expected8);
  TESTT((b * e * f) * (a * c * d * g * h), expected8);
  TESTT((b * e * g) * (a * c * d * f * h), expected8);
  TESTT((b * e * h) * (a * c * d * f * g), expected8);
  TESTT((b * f * g) * (a * c * d * e * h), expected8);
  TESTT((b * f * h) * (a * c * d * e * g), expected8);
  TESTT((b * g * h) * (a * c * d * e * f), expected8);
  TESTT((c * d * e) * (a * b * f * g * h), expected8);
  TESTT((c * d * f) * (a * b * e * g * h), expected8);
  TESTT((c * d * g) * (a * b * e * f * h), expected8);
  TESTT((c * d * h) * (a * b * e * f * g), expected8);
  TESTT((c * e * f) * (a * b * d * g * h), expected8);
  TESTT((c * e * g) * (a * b * d * f * h), expected8);
  TESTT((c * e * h) * (a * b * d * f * g), expected8);
  TESTT((c * f * g) * (a * b * d * e * h), expected8);
  TESTT((c * f * h) * (a * b * d * e * g), expected8);
  TESTT((c * g * h) * (a * b * d * e * f), expected8);
  TESTT((d * e * f) * (a * b * c * g * h), expected8);
  TESTT((d * e * g) * (a * b * c * f * h), expected8);
  TESTT((d * e * h) * (a * b * c * f * g), expected8);
  TESTT((d * f * g) * (a * b * c * e * h), expected8);
  TESTT((d * f * h) * (a * b * c * e * g), expected8);
  TESTT((d * g * h) * (a * b * c * e * f), expected8);
  TESTT((e * f * g) * (a * b * c * d * h), expected8);
  TESTT((e * f * h) * (a * b * c * d * g), expected8);
  TESTT((e * g * h) * (a * b * c * d * f), expected8);
  TESTT((f * g * h) * (a * b * c * d * e), expected8);
  TESTT((a * b * c * d) * (e * f * g * h), expected8);
  TESTT((a * b * c * e) * (d * f * g * h), expected8);
  TESTT((a * b * c * f) * (d * e * g * h), expected8);
  TESTT((a * b * c * g) * (d * e * f * h), expected8);
  TESTT((a * b * c * h) * (d * e * f * g), expected8);
  TESTT((a * b * d * e) * (c * f * g * h), expected8);
  TESTT((a * b * d * f) * (c * e * g * h), expected8);
  TESTT((a * b * d * g) * (c * e * f * h), expected8);
  TESTT((a * b * d * h) * (c * e * f * g), expected8);
  TESTT((a * b * e * f) * (c * d * g * h), expected8);
  TESTT((a * b * e * g) * (c * d * f * h), expected8);
  TESTT((a * b * e * h) * (c * d * f * g), expected8);
  TESTT((a * b * f * g) * (c * d * e * h), expected8);
  TESTT((a * b * f * h) * (c * d * e * g), expected8);
  TESTT((a * b * g * h) * (c * d * e * f), expected8);
  TESTT((a * c * d * e) * (b * f * g * h), expected8);
  TESTT((a * c * d * f) * (b * e * g * h), expected8);
  TESTT((a * c * d * g) * (b * e * f * h), expected8);
  TESTT((a * c * d * h) * (b * e * f * g), expected8);
  TESTT((a * c * e * f) * (b * d * g * h), expected8);
  TESTT((a * c * e * g) * (b * d * f * h), expected8);
  TESTT((a * c * e * h) * (b * d * f * g), expected8);
  TESTT((a * c * f * g) * (b * d * e * h), expected8);
  TESTT((a * c * f * h) * (b * d * e * g), expected8);
  TESTT((a * c * g * h) * (b * d * e * f), expected8);
  TESTT((a * d * e * f) * (b * c * g * h), expected8);
  TESTT((a * d * e * g) * (b * c * f * h), expected8);
  TESTT((a * d * e * h) * (b * c * f * g), expected8);
  TESTT((a * d * f * g) * (b * c * e * h), expected8);
  TESTT((a * d * f * h) * (b * c * e * g), expected8);
  TESTT((a * d * g * h) * (b * c * e * f), expected8);
  TESTT((a * e * f * g) * (b * c * d * h), expected8);
  TESTT((a * e * f * h) * (b * c * d * g), expected8);
  TESTT((a * e * g * h) * (b * c * d * f), expected8);
  TESTT((a * f * g * h) * (b * c * d * e), expected8);
  TESTT((b * c * d * e) * (a * f * g * h), expected8);
  TESTT((b * c * d * f) * (a * e * g * h), expected8);
  TESTT((b * c * d * g) * (a * e * f * h), expected8);
  TESTT((b * c * d * h) * (a * e * f * g), expected8);
  TESTT((b * c * e * f) * (a * d * g * h), expected8);
  TESTT((b * c * e * g) * (a * d * f * h), expected8);
  TESTT((b * c * e * h) * (a * d * f * g), expected8);
  TESTT((b * c * f * g) * (a * d * e * h), expected8);
  TESTT((b * c * f * h) * (a * d * e * g), expected8);
  TESTT((b * c * g * h) * (a * d * e * f), expected8);
  TESTT((b * d * e * f) * (a * c * g * h), expected8);
  TESTT((b * d * e * g) * (a * c * f * h), expected8);
  TESTT((b * d * e * h) * (a * c * f * g), expected8);
  TESTT((b * d * f * g) * (a * c * e * h), expected8);
  TESTT((b * d * f * h) * (a * c * e * g), expected8);
  TESTT((b * d * g * h) * (a * c * e * f), expected8);
  TESTT((b * e * f * g) * (a * c * d * h), expected8);
  TESTT((b * e * f * h) * (a * c * d * g), expected8);
  TESTT((b * e * g * h) * (a * c * d * f), expected8);
  TESTT((b * f * g * h) * (a * c * d * e), expected8);
  TESTT((c * d * e * f) * (a * b * g * h), expected8);
  TESTT((c * d * e * g) * (a * b * f * h), expected8);
  TESTT((c * d * e * h) * (a * b * f * g), expected8);
  TESTT((c * d * f * g) * (a * b * e * h), expected8);
  TESTT((c * d * f * h) * (a * b * e * g), expected8);
  TESTT((c * d * g * h) * (a * b * e * f), expected8);
  TESTT((c * e * f * g) * (a * b * d * h), expected8);
  TESTT((c * e * f * h) * (a * b * d * g), expected8);
  TESTT((c * e * g * h) * (a * b * d * f), expected8);
  TESTT((c * f * g * h) * (a * b * d * e), expected8);
  TESTT((d * e * f * g) * (a * b * c * h), expected8);
  TESTT((d * e * f * h) * (a * b * c * g), expected8);
  TESTT((d * e * g * h) * (a * b * c * f), expected8);
  TESTT((d * f * g * h) * (a * b * c * e), expected8);
  TESTT((e * f * g * h) * (a * b * c * d), expected8);
  TESTT((a * b * c * d * e) * (f * g * h), expected8);
  TESTT((a * b * c * d * f) * (e * g * h), expected8);
  TESTT((a * b * c * d * g) * (e * f * h), expected8);
  TESTT((a * b * c * d * h) * (e * f * g), expected8);
  TESTT((a * b * c * e * f) * (d * g * h), expected8);
  TESTT((a * b * c * e * g) * (d * f * h), expected8);
  TESTT((a * b * c * e * h) * (d * f * g), expected8);
  TESTT((a * b * c * f * g) * (d * e * h), expected8);
  TESTT((a * b * c * f * h) * (d * e * g), expected8);
  TESTT((a * b * c * g * h) * (d * e * f), expected8);
  TESTT((a * b * d * e * f) * (c * g * h), expected8);
  TESTT((a * b * d * e * g) * (c * f * h), expected8);
  TESTT((a * b * d * e * h) * (c * f * g), expected8);
  TESTT((a * b * d * f * g) * (c * e * h), expected8);
  TESTT((a * b * d * f * h) * (c * e * g), expected8);
  TESTT((a * b * d * g * h) * (c * e * f), expected8);
  TESTT((a * b * e * f * g) * (c * d * h), expected8);
  TESTT((a * b * e * f * h) * (c * d * g), expected8);
  TESTT((a * b * e * g * h) * (c * d * f), expected8);
  TESTT((a * b * f * g * h) * (c * d * e), expected8);
  TESTT((a * c * d * e * f) * (b * g * h), expected8);
  TESTT((a * c * d * e * g) * (b * f * h), expected8);
  TESTT((a * c * d * e * h) * (b * f * g), expected8);
  TESTT((a * c * d * f * g) * (b * e * h), expected8);
  TESTT((a * c * d * f * h) * (b * e * g), expected8);
  TESTT((a * c * d * g * h) * (b * e * f), expected8);
  TESTT((a * c * e * f * g) * (b * d * h), expected8);
  TESTT((a * c * e * f * h) * (b * d * g), expected8);
  TESTT((a * c * e * g * h) * (b * d * f), expected8);
  TESTT((a * c * f * g * h) * (b * d * e), expected8);
  TESTT((a * d * e * f * g) * (b * c * h), expected8);
  TESTT((a * d * e * f * h) * (b * c * g), expected8);
  TESTT((a * d * e * g * h) * (b * c * f), expected8);
  TESTT((a * d * f * g * h) * (b * c * e), expected8);
  TESTT((a * e * f * g * h) * (b * c * d), expected8);
  TESTT((b * c * d * e * f) * (a * g * h), expected8);
  TESTT((b * c * d * e * g) * (a * f * h), expected8);
  TESTT((b * c * d * e * h) * (a * f * g), expected8);
  TESTT((b * c * d * f * g) * (a * e * h), expected8);
  TESTT((b * c * d * f * h) * (a * e * g), expected8);
  TESTT((b * c * d * g * h) * (a * e * f), expected8);
  TESTT((b * c * e * f * g) * (a * d * h), expected8);
  TESTT((b * c * e * f * h) * (a * d * g), expected8);
  TESTT((b * c * e * g * h) * (a * d * f), expected8);
  TESTT((b * c * f * g * h) * (a * d * e), expected8);
  TESTT((b * d * e * f * g) * (a * c * h), expected8);
  TESTT((b * d * e * f * h) * (a * c * g), expected8);
  TESTT((b * d * e * g * h) * (a * c * f), expected8);
  TESTT((b * d * f * g * h) * (a * c * e), expected8);
  TESTT((b * e * f * g * h) * (a * c * d), expected8);
  TESTT((c * d * e * f * g) * (a * b * h), expected8);
  TESTT((c * d * e * f * h) * (a * b * g), expected8);
  TESTT((c * d * e * g * h) * (a * b * f), expected8);
  TESTT((c * d * f * g * h) * (a * b * e), expected8);
  TESTT((c * e * f * g * h) * (a * b * d), expected8);
  TESTT((d * e * f * g * h) * (a * b * c), expected8);
  TESTT((a * b * c * d * e * f) * (g * h), expected8);
  TESTT((a * b * c * d * e * g) * (f * h), expected8);
  TESTT((a * b * c * d * e * h) * (f * g), expected8);
  TESTT((a * b * c * d * f * g) * (e * h), expected8);
  TESTT((a * b * c * d * f * h) * (e * g), expected8);
  TESTT((a * b * c * d * g * h) * (e * f), expected8);
  TESTT((a * b * c * e * f * g) * (d * h), expected8);
  TESTT((a * b * c * e * f * h) * (d * g), expected8);
  TESTT((a * b * c * e * g * h) * (d * f), expected8);
  TESTT((a * b * c * f * g * h) * (d * e), expected8);
  TESTT((a * b * d * e * f * g) * (c * h), expected8);
  TESTT((a * b * d * e * f * h) * (c * g), expected8);
  TESTT((a * b * d * e * g * h) * (c * f), expected8);
  TESTT((a * b * d * f * g * h) * (c * e), expected8);
  TESTT((a * b * e * f * g * h) * (c * d), expected8);
  TESTT((a * c * d * e * f * g) * (b * h), expected8);
  TESTT((a * c * d * e * f * h) * (b * g), expected8);
  TESTT((a * c * d * e * g * h) * (b * f), expected8);
  TESTT((a * c * d * f * g * h) * (b * e), expected8);
  TESTT((a * c * e * f * g * h) * (b * d), expected8);
  TESTT((a * d * e * f * g * h) * (b * c), expected8);
  TESTT((b * c * d * e * f * g) * (a * h), expected8);
  TESTT((b * c * d * e * f * h) * (a * g), expected8);
  TESTT((b * c * d * e * g * h) * (a * f), expected8);
  TESTT((b * c * d * f * g * h) * (a * e), expected8);
  TESTT((b * c * e * f * g * h) * (a * d), expected8);
  TESTT((b * d * e * f * g * h) * (a * c), expected8);
  TESTT((c * d * e * f * g * h) * (a * b), expected8);
  TESTT((a * b * c * d * e * f * g) * h, expected8);
  TESTT((a * b * c * d * e * f * h) * g, expected8);
  TESTT((a * b * c * d * e * g * h) * f, expected8);
  TESTT((a * b * c * d * f * g * h) * e, expected8);
  TESTT((a * b * c * e * f * g * h) * d, expected8);
  TESTT((a * b * d * e * f * g * h) * c, expected8);
  TESTT((a * c * d * e * f * g * h) * b, expected8);
  TESTT((b * c * d * e * f * g * h) * a, expected8);

  // Test exponentiation.
  TESTS((a * (b^-two) * c * d) * ((b^-two) * d), "a * b^-4 * c * d^2");
  TESTS((a * (b^-two) * c * d) * ((b^-one) * d), "a * b^-3 * c * d^2");
  TESTS((a * (b^-two) * c * d) * ((b^zero) * d), "a * b^-2 * c * d^2"); // ⎫ same
  TESTS((a * (b^-two) * c * d) * ( one     * d), "a * b^-2 * c * d^2"); // ⎭
  TESTS((a * (b^-two) * c * d) * ((b^one)  * d), "a * b^-1 * c * d^2"); // ⎫ same
  TESTS((a * (b^-two) * c * d) * ( b       * d), "a * b^-1 * c * d^2"); // ⎭
  TESTS((a * (b^-two) * c * d) * ((b^two)  * d), "a * c * d^2");

  TESTS((a * (b^-one) * c * d) * ((b^-two) * d), "a * b^-3 * c * d^2");
  TESTS((a * (b^-one) * c * d) * ((b^-one) * d), "a * b^-2 * c * d^2");
  TESTS((a * (b^-one) * c * d) * ((b^zero) * d), "a * b^-1 * c * d^2"); // ⎫ same
  TESTS((a * (b^-one) * c * d) * ( one     * d), "a * b^-1 * c * d^2"); // ⎭
  TESTS((a * (b^-one) * c * d) * ((b^one)  * d), "a * c * d^2");        // ⎫ same
  TESTS((a * (b^-one) * c * d) * ( b       * d), "a * c * d^2");        // ⎭
  TESTS((a * (b^-one) * c * d) * ((b^two)  * d), "a * b * c * d^2");

  TESTS((a * b * c * d) * ((b^-two) * d), "a * b^-1 * c * d^2");
  TESTS((a * b * c * d) * ((b^-one) * d), "a * c * d^2");
  TESTS((a * b * c * d) * ((b^zero) * d), "a * b * c * d^2");           // ⎫ same
  TESTS((a * b * c * d) * ( one     * d), "a * b * c * d^2");           // ⎭
  TESTS((a * b * c * d) * ((b^one)  * d), "a * b^2 * c * d^2");         // ⎫ same
  TESTS((a * b * c * d) * ( b       * d), "a * b^2 * c * d^2");         // ⎭
  TESTS((a * b * c * d) * ((b^two)  * d), "a * b^3 * c * d^2");

  TESTS((a * (b^two) * c * d) * ((b^-two) * d), "a * c * d^2");
  TESTS((a * (b^two) * c * d) * ((b^-one) * d), "a * b * c * d^2");
  TESTS((a * (b^two) * c * d) * ((b^zero) * d), "a * b^2 * c * d^2");   // ⎫ same
  TESTS((a * (b^two) * c * d) * ( one     * d), "a * b^2 * c * d^2");   // ⎭
  TESTS((a * (b^two) * c * d) * ((b^one)  * d), "a * b^3 * c * d^2");   // ⎫ same
  TESTS((a * (b^two) * c * d) * ( b       * d), "a * b^3 * c * d^2");   // ⎭
  TESTS((a * (b^two) * c * d) * ((b^two)  * d), "a * b^4 * c * d^2");

  TESTS((a * b) * ((a^-two) * (c^-one)), "a^-1 * b * c^-1");

  // Test combining constants.
  TESTS(constant(42) * (constant(1, 42) * a), "a");
  TESTS(constant(42) * (constant(0, 42) * a), "0");
  TESTS(constant(42) * (constant(2, 42) * a), "2 * a");
  TESTS((a * constant(42) * c * d) * (constant(1, 42) * d), "a * c * d^2");
  TESTS((a * constant(42) * c * d) * (constant(0, 42) * d), "0");
  TESTS((a * constant(42) * c * d) * (constant(2, 42) * d), "2 * a * c * d^2");
  TESTS((two * a) * zero, "0");
  TESTS((x * y) * (x^-one), "y");

  // Test division.
  TESTS(a / b, "a * b^-1");
  TESTS((a * b) / (b * c), "a * c^-1");
  TESTS((a * b) / (a * c), "b * c^-1");
  TESTS((a * b) / ((b^two) * c), "a * b^-1 * c^-1");
  TESTS(one / ((a^two) * c), "a^-2 * c^-1");
  TESTS((a * constant(1, 2) * b) / -((a^two) * c), "-1/2 * a^-1 * b * c^-1");
  TESTS((z^two) * z, "z^3");
  TESTS(x * z / y * (x^two) / (z^(-two)), "x^3 * y^-1 * z^3");
  TESTS(x * z / y * (x^two) / -(z^(-two)), "-(x^3 * y^-1 * z^3)");
  TESTS(x * two * z / y * (x^two) / -(z^(-two)), "-2 * x^3 * y^-1 * z^3");
  TESTS(y * x * two * z / y * (x^two) / -(z^(-two)) / z, "-2 * x^3 * z^2");

  // Test printing of parenthesis.
  auto& x_plus_minus_y = x + -y;
  TESTS(x_plus_minus_y, "x - y");
  auto& x_squared = x^2;
  auto& negation_of_power = -x_squared;
  TESTS(negation_of_power, "-(x^2)");
  TESTS(x * (y^two), "x * y^2");
  TESTS(x + (y^two), "x + y^2");
  TESTS(minus__x_times_y, "-(x * y)");
  TESTS(x * y * z, "x * y * z");
  TESTS(w + x * y, "w + x * y");
  TESTS(-(x + y), "-x - y");
  TESTS(w * (x + y), "w * x + w * y");
  TESTS(w + x + y, "w + x + y");
  TESTS((x + y) * w, "w * x + w * y");

  // Addition of one symbol.
  TESTS(zero + d, "d");
  TESTS(a + zero, "a");
  TESTS(a + a, "2 * a");
  TESTS(a + two * a, "3 * a");
  TESTS(two * a + two * a, "4 * a");
  TESTS(two * a + (-one) * a, "a");
  TESTS(a + (-one) * a, "0");
  TESTS((-one) * c + c + d, "d");
  TESTS(d + (-one) * c, "-c + d");
  TESTS((-c + d) + c, "d");

  // Addition of two symbols.
  TESTS(a + b, "a + b");
  TESTS(b + a, "a + b");
  TESTS(a - b, "a - b");
  TESTS(a + two * b - c + two * (c^two) + (-one) * b - a, "b - c + 2 * c^2");

  // Addition of three symbols.
  TESTS((a + b) + (c + d), "a + b + c + d");

  auto& some_constant = Constant::realize(3);
  auto& some_symbol = y;
  auto& some_power = x^2;
  auto& some_low_product = some_constant * some_symbol;
  auto& some_product = some_symbol * z;
  auto & some_high_product = some_power * some_product;

  // Test Sum::is_less.
  // Compare constants.
  ASSERT(!Sum::is_less(constant(3, 1), constant(10, 3))); // All constants must be treated as equal.
  ASSERT(!Sum::is_less(constant(10, 3), constant(3, 1))); // All constants must be treated as equal.
  ASSERT(!Sum::is_less(some_constant, some_constant)); // A constant is not less than itself.
  // Compare symbols.
  ASSERT(Sum::is_less(constant(1000, 1), x)); // Any constant should always be considered less than a symbol.
  ASSERT(!Sum::is_less(x, constant(-1000, 1))); // A symbol is never less than a constant.
  ASSERT(Sum::is_less(x, y)); // Symbols should compare like their Id.
  ASSERT(!Sum::is_less(y, x)); // Symbols should compare like their Id.
  ASSERT(!Sum::is_less(y, y)); // A symbol is not less than itself.
  // Compare powers.
  ASSERT(Sum::is_less(some_constant, some_power)); // Any constant should always be considered less than a power.
  ASSERT(!Sum::is_less(some_power, some_constant)); // Any power is never less than a constant.
  ASSERT(Sum::is_less(some_symbol, some_power)); // Any symbol should always be considered less than a power.
  ASSERT(!Sum::is_less(some_power, some_symbol)); // Any power is never less than a symbol.
  ASSERT(Sum::is_less(x^3, x^constant(10, 3))); // Powers of equal symbol should be compared by their exponent.
  ASSERT(!Sum::is_less(x^constant(10, 3), x^3)); // Powers of equal symbol should be compared by their exponent.
  ASSERT(Sum::is_less(x^3, y^3)); // Powers should first be compared by base symbol.
  ASSERT(!Sum::is_less(y^3, x^3)); // Powers should first be compared by base symbol.
  ASSERT(!Sum::is_less(some_power, some_power)); // A power is not less than itself.
  // Compare products.
  ASSERT(Sum::is_less(some_constant, some_low_product)); // Constants are always less.
  ASSERT(!Sum::is_less(some_low_product, some_constant)); // Constants are always less.
  ASSERT(!Sum::is_less(some_symbol, some_low_product)); // some_low_product is a constant times some_symbol.
  ASSERT(!Sum::is_less(some_low_product, some_symbol)); // some_low_product is a constant times some_symbol.
  ASSERT(!Sum::is_less(some_power, some_low_product)); // some_low_product is a constant times some_symbol.
  ASSERT(Sum::is_less(some_low_product, some_power)); // some_low_product is a constant times some_symbol.
  ASSERT(Sum::is_less(some_low_product, some_high_product)); //
  ASSERT(!Sum::is_less(some_high_product, some_low_product)); //

  // Lets define four types that compare like K < L < M < N.
  auto& k = some_constant;
  auto& l = x;
  auto& m = y^2;
  auto& n = z^3;

  // Any Product<X, Y> always must have X < Y; therefore we have the following possibilities:
  auto& KL = k * l;
  auto& KM = k * m;
  auto& KN = k * n;
  auto& LM = l * m;
  auto& LN = l * n;
  auto& MN = m * n;
  // Test that we can't compile the other possibilities.
#if WONT_COMPILE
  Product<K, K> test1;  // error: static assertion [...]: The second factor of a Product can only be a Symbol, Power or another Product.
  Product<L, K> test2;  // Idem
  Product<M, K> test3;  // Idem
  Product<L, L> test4;  // error: static assertion [...]: The first argument of a Product must be less than the second argument.
  Product<M, L> test5;  // Idem
  Product<M, M> test6;  // Idem
#endif

  compare_equal1(KL);
  compare_less(KL, KM);
  compare_less(KL, KN);
  compare_less(KL, LM);
  compare_less(KL, LN);
  compare_less(KL, MN);

  compare_equal1(KM);
  compare_less(KM, KN);
  compare_less(KM, LM);
  compare_less(KM, LN);
  compare_less(KM, MN);

  compare_equal1(KN);
  compare_less(KN, LM);
  compare_less(KN, LN);
  compare_less(KN, MN);

  compare_equal1(LM);
  compare_less(LM, LN);
  compare_less(LM, MN);

  compare_equal1(LN);
  compare_less(LN, MN);

  compare_equal1(MN);

  TESTS((-a + two * b - three_halfs * c) * (x - two * y + three_halfs * c - (b^two)),
      "-2 * b^3 - 9/4 * c^2 - 3/2 * a * c - a * x + 2 * a * y + a * b^2 + 3 * b * c + 2 * b * x - 4 * b * y - 3/2 * c * x + 3 * c * y + 3/2 * b^2 * c");
  TESTS(x^(one + two), "x^3");

  TESTS((-b + (((b^two) - constant(4, 1) * a * c)^constant(1, 2))) / (two * a), "-1/2 * a^-1 * b + 1/2 * a^-1 * (b^2 - 4 * a * c)^(1/2)");
  TESTS((a + b + c + d) * (x + y + z), "a * x + a * y + a * z + b * x + b * y + b * z + c * x + c * y + c * z + d * x + d * y + d * z");
  TESTS(((a + b)^constant(3)) * x, "x * (a + b)^3");

  // Test Exponentiation.
  TESTS((constant(3, 7)^two), "9/49");
  TESTS((constant(3, 7)^constant(3)), "27/343");
  TESTS((constant(3, 7)^constant(5)), "243/16807");
  TESTS(x^two, "x^2");
  TESTS((x^three_halfs)^two, "x^3");
  TESTS(((x^constant(3, 1))^constant(1, 3)), "x");
  TESTS((((a + b)^constant(3, 1))^constant(1, 3)), "a + b");
  TESTS((constant(42) * x * y)^two, "1764 * x^2 * y^2");
  TESTS((x * (y^three_halfs))^two, "x^2 * y^3");
  TESTS((a + b)^two, "(a + b)^2");
  TESTS(((a + b)^two)^three_halfs, "(a + b)^3");
  TESTS(((a + two * b)^three_halfs) / ((x - (y^two))^two), "(a + 2 * b)^(3/2) * (x - y^2)^-2");

  // Test Division.
  TESTS(x / (y^constant(-1)), "x * y");
  TESTS(x / ((a + b)^constant(-1)), "a * x + b * x");

  auto& formula1 = (x + y)^three_halfs;
  TESTS(formula1, "(x + y)^(3/2)");
  auto& formula2 = constant(-1, 2) * formula1;
  TESTS(formula2, "-1/2 * (x + y)^(3/2)");
  TESTS(((x + y)^three_halfs) * (constant(7) * x - constant(1, 7) * y), "7 * x * (x + y)^(3/2) - 1/7 * y * (x + y)^(3/2)");

  TESTS((-two * x) * (constant(5) * ((y - (x^two))^constant(4))), "-10 * x * (y - x^2)^4");

  // Test Sin.
  TESTS(constant(2) * sin(x), "2 * sin(x)");

  // Testing is_less_Sum.
  Constant const& constant_m1 = constant(37, 41); // 0.902439024...
  Constant const& constant2 = constant(115, 127); // 0.905511811...

  // Two constants compare equal.
  compare_equal_ConstantsCompareEqual(constant_m1, constant2);
  // Unless they don't.
  compare_less_ConstantsCompareNotEqual(constant_m1, constant2);

  Expression const& symbol0 = a;
  Expression const& symbol1 = b;

  // Two symbols compare by Id.
  compare_less(symbol0, symbol1);

  Expression const& power0 = symbol0 ^ constant_m1;
  Expression const& product_c0 = constant_m1 * symbol0;
  Expression const& sin0 = sin(symbol0);
  Expression const& multiplication_cS0 = constant_m1 * sin0;
  Expression const& cos0 = cos(symbol0);
  Expression const& sum_c0 = constant_m1 + symbol0;
  Expression const& sum_c1 = constant_m1 + symbol1;
  Expression const& exponentiation = sum_c0 ^ constant_m1;

  // A constant compares less than anything else.
  compare_less(constant_m1, symbol0);
  compare_less(constant_m1, power0);
  compare_less(constant_m1, product_c0);
  compare_less(constant_m1, exponentiation);
  compare_less(constant_m1, multiplication_cS0);
  compare_less(constant_m1, sin0);
  compare_less(constant_m1, cos0);
  compare_less(constant_m1, sum_c0);

  Expression const& product_c1 = constant_m1 * symbol1;
  Expression const& product_01 = symbol0 * symbol1;
  Expression const& sin1 = sin(symbol1);
  Expression const& multiplication_0E = symbol0 * exponentiation;

  // A symbol compares less than anything else (except constants, etc - which is from now on considered trivial)
  // except other symbols that have a smaller Id (including a Product/Multiplication of a constant with such symbol).
  compare_equal(symbol0, symbol0);
  compare_less(symbol0, symbol1);
  compare_less(symbol0, power0);
  compare_equal_ConstantsCompareEqual(symbol0, product_c0);
  compare_less_ConstantsCompareNotEqual(symbol0, product_c0);
  compare_less(symbol0, product_c1);
  compare_less(symbol0, product_01);
  compare_less(symbol0, exponentiation);
  compare_less(symbol0, multiplication_cS0);
  compare_less(symbol0, multiplication_0E);
  compare_less(symbol0, sin0);
  compare_less(symbol0, cos0);
  compare_less(symbol0, sum_c0);

  Expression const& power0b = symbol0 ^ constant2;
  Expression const& power1 = symbol1 ^ constant_m1;
  Expression const& product_d0 = Product::multiply(constant2, power0);
  Expression const& product_d0b = Product::multiply(constant2, power0b);
  Expression const& product_d1 = Product::multiply(constant2, power1);

  // A power compares less than anything "larger than" a Power.
  // It compares less than another Power if either the id of the base symbol is less, or
  // the base is the same and the exponent is less.
  compare_equal(power0, power0);
  compare_less(power0, power0b);
  compare_less(power0b, power1);
  compare_equal_ConstantsCompareEqual(power0, product_d0);
  compare_less_ConstantsCompareNotEqual(power0, product_d0);
  compare_less(power0, product_d0b);
  compare_less(power0b, product_d1);
  compare_less(power0, product_01);
  compare_less(power0, exponentiation);
  compare_less(power0, multiplication_cS0);
  compare_less(power0, multiplication_0E);
  compare_less(power0, sin0);
  compare_less(power0, cos0);
  compare_less(power0, sum_c0);

  // A product compares less than anything "larger than" a Product.
  // It compares less than another Product if the first argument compares less than the first argument of the second Product,
  // or when those are equal and the second argument compares less than the second argument of the second Product.
  compare_equal(product_01, product_01);
  compare_equal_ConstantsCompareEqual(Product::multiply(constant_m1, product_01), product_01);
  compare_less_ConstantsCompareNotEqual(Product::multiply(constant_m1, product_01), product_01);
  compare_less(Product::multiply(constant_m1, exponentiation), product_01);
  compare_less(product_01, multiplication_0E);
  compare_less_ConstantsCompareEqual(product_01, Product::multiply(constant_m1, multiplication_0E));
  compare_less_ConstantsCompareNotEqual(Product::multiply(constant_m1, multiplication_0E), product_01);
  compare_less(product_01, sin0);
  compare_less_ConstantsCompareEqual(product_01, Product::multiply(constant_m1, sin0));
  compare_less_ConstantsCompareNotEqual(Product::multiply(constant_m1, sin0), product_01);
  compare_less(product_01, cos0);
  compare_less_ConstantsCompareEqual(product_01, Product::multiply(constant_m1, cos0));
  compare_less_ConstantsCompareNotEqual(Product::multiply(constant_m1, cos0), product_01);
  compare_less(product_01, sum_c0);
  compare_less(product_01, Product::multiply(constant_m1, sum_c0));

  compare_equal_ConstantsCompareEqual(product_c0, symbol0);
  compare_less_ConstantsCompareNotEqual(symbol0, product_c0);
  compare_less_ConstantsCompareEqual(product_c0, symbol1);
  compare_less_ConstantsCompareNotEqual(symbol1, product_c0);
  compare_less_ConstantsCompareEqual(product_c0, power0);
  compare_less_ConstantsCompareNotEqual(power0, product_c0);
  compare_equal(product_c0, product_c0);
  compare_less(product_c0, product_c1);
  compare_less(product_c0, product_01);
  compare_less_ConstantsCompareEqual(product_c0, exponentiation);
  compare_less_ConstantsCompareNotEqual(exponentiation, product_c0);
  compare_less(product_c0, multiplication_cS0);
  compare_less(product_c0, multiplication_0E);
  compare_less(product_c0, sin0);
  compare_less(product_c0, cos0);
  compare_less(product_c0, sum_c0);

  compare_equal_ConstantsCompareEqual(Product::multiply(constant_m1, power0), power0);
  compare_less_ConstantsCompareNotEqual(power0, Product::multiply(constant_m1, power0));
  compare_less_ConstantsCompareEqual(Product::multiply(constant_m1, power0), power0b);
  compare_less_ConstantsCompareNotEqual(power0b, Product::multiply(constant_m1, power0));
  compare_less_ConstantsCompareEqual(Product::multiply(constant_m1, power0b), power1);
  compare_less_ConstantsCompareNotEqual(power1, Product::multiply(constant_m1, power0b));
  compare_equal_ConstantsCompareEqual(Product::multiply(constant_m1, power0), product_d0);
  compare_less_ConstantsCompareNotEqual(Product::multiply(constant_m1, power0), product_d0);
  compare_less(Product::multiply(constant_m1, power0), product_d0b);
  compare_less(Product::multiply(constant_m1, power0b), product_d1);
  compare_less(Product::multiply(constant_m1, power0), product_01);
  compare_less_ConstantsCompareEqual(Product::multiply(constant_m1, power0), exponentiation);
  compare_less_ConstantsCompareNotEqual(exponentiation, Product::multiply(constant_m1, power0));
  compare_less(Product::multiply(constant_m1, power0), multiplication_cS0);
  compare_less(Product::multiply(constant_m1, power0), multiplication_0E);
  compare_less(Product::multiply(constant_m1, power0), sin0);
  compare_less(Product::multiply(constant_m1, power0), cos0);
  compare_less(Product::multiply(constant_m1, power0), sum_c0);

  Expression const& multiplication_cS1 = Product::multiply(constant_m1, sin1);

  //Expression const& exponentiation = Power::make_power(sum_c0, constant_m1);
  Expression const& exponentiation0b = Power::make_power(sum_c0, constant2);
  Expression const& exponentiation1 = Power::make_power(sum_c1, constant_m1);

  // An Exponentiation compares less than anything "larger than" an Exponentiation.
  // It compares less than another Exponentiation if either the base is less, or
  // the base is the same and the exponent is less.
  compare_equal(exponentiation, exponentiation);
  compare_less(exponentiation, product_01);
  compare_less(exponentiation, multiplication_0E);
  compare_less(exponentiation, sin0);
  compare_less(exponentiation, cos0);
  compare_less(exponentiation, sum_c0);
  compare_less(exponentiation, exponentiation0b);
  compare_less(exponentiation0b, exponentiation1);

  compare_equal_ConstantsCompareEqual(exponentiation, Product::multiply(constant_m1, exponentiation));
  compare_less_ConstantsCompareEqual(exponentiation, Product::multiply(constant_m1, multiplication_0E));
  compare_less_ConstantsCompareEqual(exponentiation, Product::multiply(constant_m1, sin0));
  compare_less_ConstantsCompareEqual(exponentiation, Product::multiply(constant_m1, cos0));
  compare_less_ConstantsCompareEqual(exponentiation, Product::multiply(constant_m1, sum_c0));
  compare_less_ConstantsCompareEqual(exponentiation, Product::multiply(constant_m1, exponentiation0b));
  compare_less_ConstantsCompareEqual(exponentiation0b, Product::multiply(constant_m1, exponentiation1));

  // A Multiplication compares less than anything "larger than" a Multiplication.
  // It compares less than another Multiplication if either the first argument is
  // less, or the first argument is equal and the second argument is less.
  compare_equal(multiplication_cS0, multiplication_cS0);
  compare_equal(multiplication_0E, multiplication_0E);
  compare_less(multiplication_0E, sin0);
  compare_less(multiplication_0E, cos0);
  compare_less(multiplication_0E, sum_c0);
  compare_less(multiplication_cS0, multiplication_cS1);

  // Symbol < Power < get_nonconstant_factor_t<Product> < Exponentiation < get_nonconstant_factor_t<Multiplication> < Sin < Cos < Sum
  compare_equal(sin0, sin0);
  compare_less(sin0, cos0);
  compare_less(sin0, sum_c0);
  compare_equal(cos0, cos0);
  compare_less(cos0, sum_c0);

  compare_equal(sum_c0, sum_c0);

  //==========================================================================
  // Old testsuite:

  test_derivative(two, x);
  test_derivative(x, x);
  test_derivative(y, x);

  x = 13;
  y = 42;

  Dout(dc::notice, "x = " << x.evaluate());
  Dout(dc::notice, "y = " << y.evaluate());

  Symbol const& y_again = Symbol::realize("y");
  Dout(dc::notice, "y_again = " << y_again.evaluate());

  Dout(dc::notice, "two.evaluate() = " << two.evaluate());

  Expression const& p = x * y;
  Dout(dc::notice, p << " [evaluate] = " << p.evaluate());
  test_derivative(p, x);
  test_derivative(p, y);

  Expression const& s = x + y;
  Dout(dc::notice, s << " [evaluate] = " << s.evaluate());
  test_derivative(s, x);
  test_derivative(s, y);

  Dout(dc::notice, "(" << s << ") * (" << s << ") = " << (s * s));
  test_derivative(s * s, x);

  Expression const& q = x^3;
  Dout(dc::notice, q << " [evaluate] = " << q.evaluate());

  auto& x_div_y = x / y;
  Dout(dc::notice, "x / y = " << x_div_y);

  test_derivative(sin(two * (x^2) * y + three * x + z), x);
  test_derivative((sin(x)^2) + (cos(x)^2), x);

  TESTS(((x + y)^two) * (((x + y + a)^two) * ((x + y)^three)), "(a + x + y)^2 * (x + y)^5");
  //TESTS(((constant(2) * a)^constant(1, 2)) * ((constant(4) * a)^constant(3, 4)), "4 * a^(5/4)");
  TESTS(((constant(2) * a)^constant(2)) * ((constant(4) * a)^constant(3)), "256 * a^5");

  Expression const& v0_div_q1_ = constant(2) * sin(y) / sin(y - x);
  Expression const& v0x_ = v0_div_q1_ * (cos(x) * c + sin(x) * a);
  Expression const& v0y_ = v0_div_q1_ * (cos(x) * d + sin(x) * b);
  Expression const& v02_ = square(v0x_) + square(v0y_);

  Dout(dc::notice, v02_);

//  Debug(Expression::dump_database());

  Function const& fn = Function::realize("fn", (x^2) + y);
  Dout(dc::notice, fn << " [" << fn.definition() << "] = " << fn.evaluate());
  test_derivative(fn, x);
  Function const& fnd = fn.derivative(x);
  Dout(dc::notice, fnd << " [" << fnd.definition() << "] = " << fnd.evaluate());

  Expression const& te = (x + (y^2))^3;
  Dout(dc::notice, te);
}
