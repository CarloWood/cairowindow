#include "sys.h"
#include "cairowindow/symbolic/symbolic.h"
#include "utils/has_print_on.h"

using namespace symbolic;

// This class defines a print_on method.
using utils::has_print_on::operator<<;

class Cubic
{
 private:
  Symbol const* x_;
  std::array<Expression const*, 4> coefficients_;

 public:
  Cubic(Symbol const& x, Expression const& c0, Expression const& c1, Expression const& c2, Expression const& c3) :
    x_(&x), coefficients_{&c0, &c1, &c2, &c3} { }

  Expression const& operator[](int i) const
  {
    ASSERT(0 <= i && i < coefficients_.size());
    return *coefficients_[i];
  }

  Cubic& operator*=(Expression const& factor)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] = &(*coefficients_[i] * factor);
    return *this;
  }

  friend Cubic operator*(Expression const& lhs, Cubic const& rhs)
  {
    Cubic result(rhs);
    result *= lhs;
    return result;
  }

  friend Cubic operator*(Cubic const& lhs, Expression const& rhs)
  {
    Cubic result(lhs);
    result *= rhs;
    return result;
  }

  Cubic& operator/=(Expression const& factor)
  {
    for (int i = 0; i < coefficients_.size(); ++i)
      coefficients_[i] = &(*coefficients_[i] / factor);
    return *this;
  }

  friend Cubic operator/(Expression const& lhs, Cubic const& rhs)
  {
    Cubic result(rhs);
    result /= lhs;
    return result;
  }

  friend Cubic operator/(Cubic const& lhs, Expression const& rhs)
  {
    Cubic result(lhs);
    result /= rhs;
    return result;
  }

  Cubic& operator+=(Expression const& term)
  {
    coefficients_[0] = &(*coefficients_[0] + term);
    return *this;
  }

  friend Cubic operator+(Expression const& lhs, Cubic const& rhs)
  {
    Cubic result(rhs);
    result += lhs;
    return result;
  }

  friend Cubic operator+(Cubic const& lhs, Expression const& rhs)
  {
    Cubic result(lhs);
    result += rhs;
    return result;
  }

  Cubic& operator-=(Expression const& term)
  {
    coefficients_[0] = &(*coefficients_[0] - term);
    return *this;
  }

  friend Cubic operator-(Expression const& lhs, Cubic const& rhs)
  {
    Cubic result(rhs);
    result -= lhs;
    return result;
  }

  friend Cubic operator-(Cubic const& lhs, Expression const& rhs)
  {
    Cubic result(lhs);
    result -= rhs;
    return result;
  }

  Cubic derivative() const
  {
    return {*x_, *coefficients_[1], 2 * *coefficients_[2], 3 * *coefficients_[3], Constant::realize(0)};
  }

  Expression const& evaluate(Expression const& x)
  {
    return *coefficients_[0] + (*coefficients_[1] + (*coefficients_[2] + *coefficients_[3] * x) * x) * x;
  }

  Cubic transform(Expression const& x_scale, Expression const& y_scale, Expression const& x_shift, Expression const& y_shift, Symbol const& xp);

  void print_on(std::ostream& os) const;
};

void Cubic::print_on(std::ostream& os) const
{
  bool saw_term = false;
  for (int i = 0; i < coefficients_.size(); ++i)
  {
    if (!Constant::is_zero(*coefficients_[i]))
    {
      if (saw_term)
        os << " + ";
      if (i == 0 || !Constant::is_one(*coefficients_[i]))
      {
        bool need_parens = needs_parens(coefficients_[i]->precedence(), Precedence::product, before);
        if (need_parens)
          os << '(';
        coefficients_[i]->print_on(os);
        if (need_parens)
          os << ')';
        if (i > 0)
          os << " * ";
      }
      if (i > 0)
        os << *x_;
      if (i == 2)
        os << "²";
      else if (i == 3)
        os << "³";
      saw_term = true;
    }
  }
  if (!saw_term)
    os << "0";
}

Cubic Cubic::transform(Expression const& x_scale, Expression const& y_scale, Expression const& x_shift, Expression const& y_shift, Symbol const& xp)
{
  DoutEntering(dc::notice, "Cubic::transform(" << x_scale << ", " << y_scale << ", " << x_shift << ", " << y_shift << ", " << xp << ")");

  // First apply the scaling:
  Expression const& c0_scaled = *coefficients_[0] * y_scale;
  Expression const& c1_scaled = *coefficients_[1] * y_scale / x_scale;
  Expression const& c2_scaled = *coefficients_[2] * y_scale / (x_scale * x_scale);
  Expression const& c3_scaled = *coefficients_[3] * y_scale / (x_scale * x_scale * x_scale);
  // Then the offset.
  Expression const& new_c0 = c0_scaled + y_shift + (-c1_scaled + (c2_scaled - c3_scaled * x_shift) * x_shift) * x_shift;
  Expression const& new_c1 = c1_scaled + (-2 * c2_scaled + 3 * c3_scaled * x_shift) * x_shift;
  Expression const& new_c2 = c2_scaled - 3 * c3_scaled * x_shift;

  return {xp, new_c0, new_c1, new_c2, c3_scaled};
}

int main()
{
  Debug(NAMESPACE_DEBUG::init());

  Constant const& zero = Constant::realize(0);
  Constant const& one = Constant::realize(1);

  // syms c0 c1 c2 x xp Ix Ex;
  Symbol const& c0 = Symbol::realize("c0");
  Symbol const& c1 = Symbol::realize("c1");
  Symbol const& c2 = Symbol::realize("c2");
  Symbol const& c3 = Symbol::realize("c3");
  Symbol const& x = Symbol::realize("x");
  Symbol const& xp = Symbol::realize("xp");
  Symbol const& xpp = Symbol::realize("xpp");
  Symbol const& scale_x = Symbol::realize("scale_x");
  Symbol const& scale_y = Symbol::realize("scale_y");
  Symbol const& shift_x = Symbol::realize("shift_x");
  Symbol const& shift_y = Symbol::realize("shift_y");

  // Monic cubic Polynomial in x:
  Cubic Px0(x, c0, c1, c2, c3);
  Dout(dc::notice, "P(x) = " << Px0);

  // Divide by the coefficient of x^3.
  Dout(dc::notice, "Divide P(x) by " << c3);
  auto Px1 = Px0.transform(one, 1 / c3, zero, zero, x);
  Dout(dc::notice, "P(x) = " << Px1);

  // First derivative:
  auto dPx1 = Px1.derivative();
  Dout(dc::notice, "P'(x) = " << dPx1);

  // Second derivative:
  auto ddPx1 = dPx1.derivative();
  Dout(dc::notice, "P\"(x) = " << ddPx1);

  Expression const& Ix1 = -c2 / (3 * c3);
  Dout(dc::notice, "Ix = " << Ix1);
  Dout(dc::notice, "Proof: P\"(Ix) = " << ddPx1.evaluate(Ix1));
  ASSERT(Constant::is_zero(ddPx1.evaluate(Ix1)));

  // Move the inflection point to x=0.
  Dout(dc::notice, "Shift P(x) Ix to the left");
  auto Px2 = Px1.transform(one, one, -Ix1, zero, x);
  Dout(dc::notice, "P(x) = " << Px2);

  auto dPx2 = Px2.derivative();
  Dout(dc::notice, "P'(x) = " << dPx2);

  Expression const& Ex2 = sqrt(((c2 / c3)^2) - 3 * c1 / c3) / 3;
  Dout(dc::notice, "Ex = " << Ex2);
  Dout(dc::notice, "Proof: P'(Ex) = " << dPx2.evaluate(Ex2));
  ASSERT(Constant::is_zero(dPx2.evaluate(Ex2)));

  Dout(dc::notice, "Scale the x-axis by a factor of Ex, so that the extrema will be at +/-1, and divide the polynomial by Ex^3 to keep it a monic polynomial.");
  auto Px3 = Px2.transform(1 / Ex2, 1 / Ex2^3, zero, zero, x);
  Dout(dc::notice, "P(x) = " << Px3);

  // Get coefficient of x.
  Expression const& C1 = Px3[1];
  Dout(dc::notice, "The coefficient of x is: " << C1);
  // This should be equal to -3.
  c1 = 6732387;
  c2 = 9138501;
  c3 = 1253698;
  Dout(dc::notice, "C1 = " << C1.evaluate());
  ASSERT(C1.evaluate() == -3);

  // Get coefficient of x^2.
  Expression const& C2 = Px3[2];
  // Is zero because Ix = 0.
  ASSERT(Constant::is_zero(C2));

  // Get coefficient of x^3.
  Expression const& C3 = Px3[3];
  ASSERT(Constant::is_one(C3));

  // Get the constant.
  Expression const& C0 = Px3[0];
  Dout(dc::notice, "C0 = " << C0);

  // This should be equal to
  //
  //            3
  //  27⋅c₀ - c₂  + 3⋅c₂⋅d
  //  ────────────────────
  //           3/2
  //          d
  //
  // where d = c₂^2 - 3 * c₁,
  // c₀ = c0/c3, c₁ = c1/c3 and c₂ = c2/c3.

  Expression const& d = ((c2 / c3)^2) - 3 * (c1 / c3);
  Expression const& C0a = (27 * (c0 / c3) - ((c2 / c3)^3) + 3 * (c2 / c3) * d) / (d^Constant::realize(3, 2));
  Dout(dc::notice, "C0a = " << C0a);

  c0 = 672938512;
  Dout(dc::notice, "C0 = " << C0.evaluate());
  Dout(dc::notice, "C0a = " << C0a.evaluate());
  ASSERT(C0.evaluate() == C0a.evaluate());

  Cubic Px4(x, C0a, Constant::realize(-3), zero, one);
  Dout(dc::notice, "which simplifies to:");
  Dout(dc::notice, "P(x) = " << Px4);

#if 0
  auto dPx3= Px3.derivative();
  Dout(dc::notice, "P'(x) = " << dPx3);
  Dout(dc::notice, "P'(1) = " << dPx3.evaluate(one));
#endif
}
