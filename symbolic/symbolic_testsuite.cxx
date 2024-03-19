#include "sys.h"
#include <sstream>
#include "Constant.h"
//#include "Negation.h"
#include "Symbol.h"
#include "Product.h"
#include "Power.h"
#include "Sum.h"
#include "debug.h"
#include <cassert>

#define TESTS(x, y) TEST(#x, x, y)

template<typename T>
void TEST(char const* what, T const& expression, std::string expect_string)
{
#ifdef SYMBOLIC_PRINTING
  std::ostringstream result;
  result << expression;
#ifdef CWDEBUG
  std::ostringstream type_os;
  type_os << NAMESPACE_DEBUG::type_name_of<T>();
  std::string type_str = type_os.str();
  size_t pos = 0;
  while ((pos = type_str.find("symbolic::", pos)) != std::string::npos)
    type_str.erase(pos, 10);
  pos = 0;
  while ((pos = type_str.find("> >", pos)) != std::string::npos)
    type_str.erase(pos + 1, 1);
  Dout(dc::notice, what << " = " << type_str << " = " << result.str());
  ASSERT(result.str() == expect_string);
#else
  std::cout << what << " = " << result.str() << std::endl;
  assert(result.str() == expect_string);
#endif
#endif  // SYMBOLIC_PRINTING
}

template<typename T, typename E>
void TESTT(T const& expression, E const& expect)
{
  static_assert(std::is_same_v<T, E>, "Failure!");
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  using namespace symbolic;

  // Constant creation.
  TEST("42 / 1", constant<42, 1>(), "42");
  TEST("-42 / 1", constant<-42, 1>(), "-42");
  TEST("42 / -1", constant<42, -1>(), "-42");
  TEST("369 / 3", constant<369, 3>(), "123");
  TEST("-369 / -3", constant<-369, -3>(), "123");
  TEST("(2 * 3 * 7 * 13 / (2 * 7 * 11 * 17)", constant<2 * 3 * 7 * 13, 2 * 7 * 11 * 17>(), "39/187");

  // Constant negation.
  TEST("-(13/42)", -constant<13, 42>(), "-13/42");
  TEST("-(-13/42)", -constant<-13, 42>(), "13/42");

  constexpr auto A = constant<2, 3 * 7>();
  constexpr auto B = constant<-5, 7 * 11>();
  constexpr auto C = constant<5, 7 * 11>();

  // Addition of two Constant's.
  TEST("A(2/21) + B(-5/77)", A + B, "1/33");

  // Subtraction of two Constant's.
  TEST("A(2/21) - B(-5/77)", A - B, "37/231");
  TEST("A(2/21) - C(5/77)", A - C, "1/33");

  // Multiplication of two Constant's.
  TEST("A(2/21) * B(-5/77)", A * B, "-10/1617");

  // Division of two Constant's.
  TEST("A(2/21) / B(-5/77)", A / B, "-22/15");

  constexpr auto C2 = constant<2, 1>();
  constexpr auto C29 = constant<29, 1>();

  // Exponentiation of a Constant.
  TEST("B(-5/77) ^ 3", B^constant<3, 1>(), "-125/456533");
  TEST("C2(2) ^ C29(29)", C2^C29, "536870912");

  // Negation of a Constant.
  TEST("-B(-5/77)", -B, "5/77");

  static constexpr auto v = make_symbol("v");
  static constexpr auto w = make_symbol("w");
  static constexpr auto x = make_symbol("x");
  static constexpr auto y = make_symbol("y");
  static constexpr auto z = make_symbol("z");

  // Negation of a Symbol.
  TEST("-x", -x, "-x");
  // Negation of a Negation.
  TEST("-(-x)", -(-x), "x");

  // Negation tests involving products.
  using x_type = std::decay_t<decltype(x)>;
  using y_type = std::decay_t<decltype(y)>;

#if 0   // Class Negation no longer exists.
  Negation<x_type> minus_x{x};
#else
  Product<Constant<-1, 1>, x_type> minus_x{constant<-1>(), x};
#endif
  TESTT(-x, minus_x);
  TESTT(-minus_x, x);

  Product<x_type, y_type> x_times_y{x, y};
#if 0   // Class Negation no longer exists.
  Negation<Product<x_type, y_type>> minus__x_times_y{x_times_y};
  Negation<y_type> minus_y{y};
#else
  Product<Constant<-1, 1>, Product<x_type, y_type>> minus__x_times_y{constant<-1>(), x_times_y};
  Product<Constant<-1, 1>, y_type> minus_y{constant<-1>(), y};
#endif
  TESTT(x * minus_y, minus__x_times_y);

  constexpr auto zero = constant<0>();
  constexpr auto one = constant<1>();
  constexpr auto two = constant<2>();
  constexpr auto three_halfs = constant<3, 2>();

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

  static constexpr auto a = make_symbol("a");
  static constexpr auto b = make_symbol("b");
  static constexpr auto c = make_symbol("c");
  static constexpr auto d = make_symbol("d");
  static constexpr auto e = make_symbol("e");
  static constexpr auto f = make_symbol("f");
  static constexpr auto g = make_symbol("g");
  static constexpr auto h = make_symbol("h");

  // Test sorting of 8 symbols during multiplication.
  auto expected8 = a * (b * (c * (d * (e * (f * (g * h))))));
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
  TESTS((a * (b^-one) * c * d) * ((b^one)  * d), "a * c * d^2");          // ⎫ same
  TESTS((a * (b^-one) * c * d) * ( b       * d), "a * c * d^2");          // ⎭
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
  TESTS(constant<42>() * (constant<1, 42>() * a), "a");
  TESTS(constant<42>() * (constant<0, 42>() * a), "0");
  TESTS(constant<42>() * (constant<2, 42>() * a), "2 * a");
  TESTS((a * constant<42>() * c * d) * (constant<1, 42>() * d), "a * c * d^2");
  TESTS((a * constant<42>() * c * d) * (constant<0, 42>() * d), "0");
  TESTS((a * constant<42>() * c * d) * (constant<2, 42>() * d), "2 * a * c * d^2");
  TESTS((two * a) * zero, "0");
  TESTS((x * y) * (x^-one), "y");

  // Test division.
  TESTS(a / b, "a * b^-1");
  TESTS((a * b) / (b * c), "a * c^-1");
  TESTS((a * b) / (a * c), "b * c^-1");
  TESTS((a * b) / ((b^two) * c), "a * b^-1 * c^-1");
  TESTS(one / ((a^two) * c), "a^-2 * c^-1");
  TESTS((a * constant<1, 2>() * b) / -((a^two) * c), "-1/2 * a^-1 * b * c^-1");
  TESTS((z^two) * z, "z^3");
  TESTS(x * z / y * (x^two) / (z^(-two)), "x^3 * y^-1 * z^3");
  TESTS(x * z / y * (x^two) / -(z^(-two)), "-(x^3 * y^-1 * z^3)");
  TESTS(x * two * z / y * (x^two) / -(z^(-two)), "-2 * x^3 * y^-1 * z^3");
  TESTS(y * x * two * z / y * (x^two) / -(z^(-two)) / z, "-2 * x^3 * z^2");

  // Test printing of parenthesis.
#if 0 // No longer allowed to create a Negation of a Negation.
  Negation<Negation<x_type>> minus_minus_x{minus_x};
  TESTS(minus_minus_x, "-(-x)");
#endif
#if 0 // Not allowed to construct the Power of a Negation.
  Power<Negation<x_type>, 2, 1> power_of_negation(minus_x);
  TESTS(power_of_negation, "(-x)^2");
#endif
#if 0   // Class Negation no longer exists.
  Product<x_type, Negation<y_type>> x_times_minus_y{x, minus_y};
  TESTS(x_times_minus_y, "x * -y");
#endif
#if 0 // Not possible to create a Ratio (which doesn't exist).
  TESTS(x_div_minus_y, "x/-y");
#endif
#if 0   // Class Negation no longer exists.
  Sum<x_type, Negation<y_type>> x_plus_minus_y{x, minus_y};
#else
  Sum<x_type, Product<Constant<-1, 1>, y_type>> x_plus_minus_y{x, minus_y};
#endif
  TESTS(x_plus_minus_y, "x + -y");
  Power<x_type, 2, 1> x_squared{x};
#if 0   // Class Negation no longer exists.
  Negation<Power<x_type, 2, 1>> negation_of_power{x_squared};
#else
  Product<Constant<-1, 1>, Power<x_type, 2, 1>> negation_of_power{constant<-1>(), x_squared};
#endif
  TESTS(negation_of_power, "-(x^2)");
#if 0 // Note allowed to create a Power of anything but a Symbol.
  Power<Power<x_type, 2, 1>, 3, 1> power_of_power{x_squared};
  TESTS(power_of_power, "(x^2)^3");
#endif
  TESTS(x * (y^two), "x * y^2");
  TESTS(x + (y^two), "x + y^2");
  TESTS(minus__x_times_y, "-(x * y)");
#if 0 // Note allowed to create a Power of anything but a Symbol.
  TESTS((x * y)^two, "(x * y)^2");
#endif
  TESTS(x * y * z, "x * y * z");
  TESTS(w + x * y, "w + x * y");
//  TESTS(-(x + y), "-(x + y)"); FIXME: uncomment once we can multiply a constant with a sum.
#if 0 // Note allowed to create a Power of anything but a Symbol.
  TESTS("(x + y)^2");
#endif
//  TESTS(w * (x + y), "w * (x + y)"); FIXME: uncomment once we can multiply a symbol with a sum.
  TESTS(w + x + y, "w + x + y");

  // Addition of one symbol.
//  TESTS(a + a, "2 * a");
//  TESTS(a + two * a, "3 * a");
//  TESTS(two * a + two * a, "4 * a");
//  TESTS(two * a + (-one) * a, "a");
//  TESTS(a + (-one) * a, "0");

  // Addition of two symbols.
  TESTS(a + b, "a + b");
  TESTS(b + a, "a + b");

  // Addition of three symbols.
  TESTS((a + b) + (c + d), "a + b + c + d");

  using some_constant_type = Constant<3, 1>;
  using some_symbol_type = y_type;
  using some_power_type = Power<x_type, 2, 1>;
  using some_low_product_type = Product<some_constant_type, some_symbol_type>;
  using some_high_product_type = Product<some_power_type, some_low_product_type>;

  // Test is_less_v.
  // Compare constants.
  static_assert(is_less_v<Constant<3, 1>, Constant<10, 3>>, "3 < 10/3 should be true.");
  static_assert(!is_less_v<Constant<10, 3>, Constant<3, 1>>, "10/3 < 3 should be false.)");
  static_assert(!is_less_v<some_constant_type, some_constant_type>, "A constant is not less than itself.");
  // Compare symbols.
  static_assert(is_less_v<Constant<1000, 1>, x_type>, "Any constant should always be considered less than a symbol.");
  static_assert(!is_less_v<x_type, Constant<-1000, 1>>, "A symbol is never less than a constant.");
  static_assert(is_less_v<x_type, y_type>, "Symbols should compare like their Id.");
  static_assert(!is_less_v<y_type, x_type>, "Symbols should compare like their Id.");
  static_assert(!is_less_v<y_type, y_type>, "A symbol is not less than itself.");
  // Compare powers.
  static_assert(is_less_v<some_constant_type, some_power_type>, "Any constant should always be considered less than a power.");
  static_assert(!is_less_v<some_power_type, some_constant_type>, "Any power is never less than a constant.");
  static_assert(is_less_v<some_symbol_type, some_power_type>, "Any symbol should always be considered less than a power.");
  static_assert(!is_less_v<some_power_type, some_symbol_type>, "Any power is never less than a symbol.");
  static_assert(is_less_v<Power<x_type, 3, 1>, Power<x_type, 10, 3>>, "Powers of equal symbol should be compared by their exponent.");
  static_assert(!is_less_v<Power<x_type, 10, 3>, Power<x_type, 3, 1>>, "Powers of equal symbol should be compared by their exponent.");
  static_assert(is_less_v<Power<x_type, 3, 1>, Power<y_type, 3, 1>>, "Powers should first be compared by base symbol.");
  static_assert(!is_less_v<Power<y_type, 3, 1>, Power<x_type, 3, 1>>, "Powers should first be compared by base symbol.");
  static_assert(!is_less_v<some_power_type, some_power_type>, "A power is not less than itself.");
  // Compare products.
  static_assert(is_less_v<some_constant_type, some_low_product_type>, "Constants are always less.");
  static_assert(!is_less_v<some_low_product_type, some_constant_type>, "Constants are always less.");
  static_assert(is_less_v<some_symbol_type, some_low_product_type>, "Symbols are always less.");
  static_assert(!is_less_v<some_low_product_type, some_symbol_type>, "Symbols are always less.");
  static_assert(is_less_v<some_power_type, some_low_product_type>, "Powers are always less.");
  static_assert(!is_less_v<some_low_product_type, some_power_type>, "Powers are always less.");
  static_assert(is_less_v<some_low_product_type, some_high_product_type>, "");
  static_assert(!is_less_v<some_high_product_type, some_low_product_type>, "");

  // Lets define three types that compare like K < L < M.
  using K = some_constant_type;
  using L = some_symbol_type;
  using M = some_power_type;

  static constexpr K k = constant<3, 1>();
  static constexpr L l = y;
  static constexpr M m{x};

  // Any Product<X, Y> always must have X < Y; therefore we have the following possibilities:
  using KL = Product<K, L>;
  using KM = Product<K, M>;
  using LM = Product<L, M>;
  // Test that we can't compile the other possibilities.
#if 0
  Product<K, K> test1{k, k};  // error: static assertion [...]: The second factor of a Product can only be a Symbol, Power or another Product.
  Product<L, K> test2{l, k};  // Idem
  Product<M, K> test3{m, k};  // Idem
#endif
  Product<L, L> test4{l, l};  // compiled!
  Product<M, L> test5{m, l};  // compiled!
  Product<M, M> test6{m, m};  //
  KL test7{k, l};
  KM test8{k, m};
  LM test9{l, m};
}
