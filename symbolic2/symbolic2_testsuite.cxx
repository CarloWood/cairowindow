#include "sys.h"
#include "Symbol.h"
#include "Constant.h"
#include "debug.h"

using namespace symbolic2;

void test_differentiate(Expression const& expression, Symbol const& symbol)
{
  Dout(dc::notice, "∂/∂" << symbol << " " << expression << " = " << expression.differentiate(symbol));
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  Constant const& two = Constant::realize(2);
  Symbol const& x = Symbol::realize("x");
  Symbol const& y = Symbol::realize("y");

  test_differentiate(two, x);
  test_differentiate(x, x);
  test_differentiate(y, x);

  x = 13;
  y = 42;

  Dout(dc::notice, "x = " << x.evaluate());
  Dout(dc::notice, "y = " << y.evaluate());

  Symbol const& z = Symbol::realize("y");
  Dout(dc::notice, "z = " << z.evaluate());

  Dout(dc::notice, "two.evaluate() = " << two.evaluate());

  //Debug(Expression::dump_database());
}
