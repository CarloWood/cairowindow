#include "sys.h"
#include <sstream>
#include "Constant.h"
#include "Negation.h"
#include "Symbol.h"
#include "Product.h"
#include "debug.h"

template<typename T>
void TEST(char const* what, T const& expression, std::string expect_string)
{
  std::stringstream oss;
  oss << expression;
  Dout(dc::notice, what << " = " << NAMESPACE_DEBUG::type_name_of<T>() << " = " << oss.str());
  ASSERT(oss.str() == expect_string);
}

template<int n>
struct Foo
{
  void print() const
  {
    std::cout << "Foo<" << n << ">\n";
  }
};

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

  static constexpr Symbol x("x");
  static constexpr Symbol y("y");

  // Negation of a Negation.
  TEST("-(-x)", -(-x), "x");

  TEST("x * y", (x * y), "x * y");
}
