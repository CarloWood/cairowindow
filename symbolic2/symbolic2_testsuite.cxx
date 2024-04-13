#include "sys.h"
#include "Product.h"
#include "Sum.h"
#include "Symbol.h"
#include "Constant.h"
#include "debug.h"

using namespace symbolic2;

void test_differentiate(Expression const& expression, Symbol const& symbol)
{
  Dout(dc::notice, "∂/∂" << symbol << " (" << expression << ") = " << expression.differentiate(symbol));
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

  Expression const& p = x * y;
  Dout(dc::notice, p << " [evaluate] = " << p.evaluate());
  test_differentiate(p, x);
  test_differentiate(p, y);

  Expression const& s = x + y;
  Dout(dc::notice, s << " [evaluate] = " << s.evaluate());
  test_differentiate(s, x);
  test_differentiate(s, y);

  Dout(dc::notice, "(" << s << ") * (" << s << ") = " << (s * s));
  test_differentiate(s * s, x);

  Expression const& q = x^3;
  Dout(dc::notice, q << " [evaluate] = " << q.evaluate());

  auto& minus_x = -x;
  Dout(dc::notice, "-x = " << minus_x);

  auto& x_div_y = x / y;
  Dout(dc::notice, "x / y = " << x_div_y);

  Debug(Expression::dump_database());
}
