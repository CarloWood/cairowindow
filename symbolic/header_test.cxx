#include "sys.h"

//#define TEST_ID_RANGE
//#define TEST_PRECEDENCE
//#define TEST_COUNTER
//#define TEST_CANONICAL_CONSTANT
//#define TEST_CONSTANT
//#define TEST_MULTIPLY_FWD
//#define TEST_EXPRESSION

//#define TEST_IS_CONSTANT
//#define TEST_SYMBOL
//#define TEST_POWER
//#define TEST_PRODUCT
//#define TEST_SUM
//#define TEST_EXPONENTIATION
//#define TEST_MULTIPLICATION
//#define TEST_SIN
//#define TEST_COS

//#define TEST_MULTIPLY
//#define TEST_ADD_EQUALS

//#define TEST_INVERT
//#define TEST_NEGATE

//#define TEST_ADD
//#define TEST_EVALUATE
//#define TEST_GET_BASE
//#define TEST_GET_EXPONENT
//#define TEST_GET_CONSTANT_FACTOR
//#define TEST_GET_NON_CONSTANT_FACTOR
//#define TEST_IS_SYMBOL
//#define TEST_IS_POWER
//#define TEST_IS_PRODUCT
//#define TEST_IS_EXPONENTIATION
//#define TEST_IS_MULTIPLICATION
//#define TEST_IS_SUM
//#define TEST_IS_SIN
//#define TEST_IS_COS
//#define TEST_IS_LESS_SUM
//#define TEST_IS_LESS_PRODUCT
//#define TEST_IS_LESS_MULTIPLICATION
//#define TEST_EXPRESSION_TRAITS
//#define TEST_OPERATORS
//#define TEST_DIFFERENTIATE
//#define TEST_SYMBOLIC

#ifdef TEST_ID_RANGE
#include "IdRange.h"
#endif
#ifdef TEST_PRECEDENCE
#include "precedence.h"
#endif
#ifdef TEST_COUNTER
#include "counter.h"
#endif
#ifdef TEST_CANONICAL_CONSTANT
#include "canonical_constant.h"
#endif
#ifdef TEST_CONSTANT
#include "Constant.h"
#endif
#ifdef TEST_MULTIPLY_FWD
#include "multiply_fwd.h"
#endif
#ifdef TEST_EXPRESSION
#include "Expression.h"
#endif

#ifdef TEST_IS_CONSTANT
#include "is_constant.h"
#include "Constant.h"
#endif

#ifdef TEST_SYMBOL
#include "Symbol.h"
#endif

#ifdef TEST_POWER
#include "Power.h"
#include "Symbol.h"     // Required for the Base::s_id of Power.
#endif
#ifdef TEST_PRODUCT
#include "Product.h"
#include "Symbol.h"
#include "Power.h"
#endif
#ifdef TEST_SUM
#include "Sum.h"
#endif
#ifdef TEST_EXPONENTIATION
#include "Exponentiation.h"
#endif
#ifdef TEST_MULTIPLICATION
#include "Multiplication.h"
#endif
#if defined(TEST_SIN) || defined(TEST_COS)
#include "functions.h"
#endif

#ifdef TEST_INVERT
#include "invert.h"
#include "Symbol.h"
#include "Power.h"
#endif
#ifdef TEST_NEGATE
#include "negate.h"
#include "Symbol.h"
#include "multiply.h"
#endif
#ifdef TEST_MULTIPLY
#include "multiply.h"
#include "Symbol.h"
#include "Power.h"
#endif
#ifdef TEST_ADD_EQUALS
#include "add_equals.h"
#include "Symbol.h"
#endif
#ifdef TEST_ADD
#include "add.h"
#include "Symbol.h"
#endif
#ifdef TEST_EVALUATE
#include "evaluate.h"
#endif
#ifdef TEST_EXPONENTIATE
#include "exponentiate.h"
#endif
#ifdef TEST_GET_BASE
#include "get_base.h"
#include "exponentiate.h"
#endif
#ifdef TEST_GET_EXPONENT
#include "get_exponent.h"
#include "exponentiate.h"
#endif
#ifdef TEST_GET_CONSTANT_FACTOR
#include "get_constant_factor.h"
#include "Constant.h"
#include "Symbol.h"
#include "Product.h"
#endif
#ifdef TEST_GET_NON_CONSTANT_FACTOR
#include "get_nonconstant_factor.h"
#include "Constant.h"
#include "Symbol.h"
#include "Power.h"
#include "Product.h"
#endif
#ifdef TEST_IS_SYMBOL
#include "is_symbol.h"
#endif
#ifdef TEST_IS_POWER
#include "is_power.h"
#endif
#ifdef TEST_IS_PRODUCT
#include "is_product.h"
#include "Symbol.h"
#include "Power.h"
#endif
#ifdef TEST_IS_EXPONENTIATION
#include "is_exponentiation.h"
#endif
#ifdef TEST_IS_MULTIPLICATION
#include "is_multiplication.h"
#endif
#ifdef TEST_IS_SUM
#include "is_sum.h"
#include "Symbol.h"
#endif
#ifdef TEST_IS_SIN
#include "is_sin.h"
#include "Symbol.h"
#endif
#ifdef TEST_IS_COS
#include "is_cos.h"
#include "Symbol.h"
#endif

#ifdef TEST_IS_LESS_SUM
#include "is_less_Sum.h"
#endif
#ifdef TEST_IS_LESS_PRODUCT
#include "is_less_Product.h"
#endif
#ifdef TEST_IS_LESS_MULTIPLICATION
#include "is_less_Multiplication.h"
#endif

#ifdef TEST_EXPRESSION_TRAITS
#include "expression_traits.h"
#endif
#ifdef TEST_OPERATORS
#include "operators.h"
#endif
#ifdef TEST_DIFFERENTIATE
#include "differentiate.h"
#endif
#ifdef TEST_SYMBOLIC
#include "symbolic.h"
#endif

#if defined(TEST_SUM) || defined(TEST_EXPONENTIATION) || defined(TEST_MULTIPLICATION) || defined(TEST_SIN) || defined(TEST_COS) || defined(TEST_EVALUATE) || defined(TEST_EXPONENTIATE) || defined(TEST_GET_BASE) || defined(TEST_GET_EXPONENT) || defined(TEST_IS_EXPONENTIATION) || defined(TEST_IS_MULTIPLICATION)
#include "Symbol.h"
#include "Power.h"
#include "Product.h"
#include "Sum.h"
#include "Multiplication.h"
#include "Exponentiation.h"
#endif

int main()
{
#ifndef TEST_COUNTER
  using namespace symbolic;
#endif

#ifdef TEST_ID_RANGE
  IdRange<13, 42> id_range;
#endif

#ifdef TEST_PRECEDENCE
  precedence test = precedence::difference;
#endif

#ifdef TEST_COUNTER
  metahack::counter<42> counter;
#endif

#ifdef TEST_CANONICAL_CONSTANT
  canonical_constant<42, 13> constant;
#endif

#ifdef TEST_CONSTANT
  auto c = constant<42, 13>();
#endif

#ifdef TEST_MULTIPLY_FWD
  // nothing.
#endif

#ifdef TEST_EXPRESSION
  // nothing.
#endif

#ifdef TEST_IS_CONSTANT
  static_assert(is_constant_v<Constant<42, 23>>, "fail");
#endif

#ifdef TEST_SYMBOL
  auto symbol = make_symbol("x");
#endif

#if defined(TEST_POWER) || defined(TEST_PRODUCT) || defined(TEST_SUM) || defined(TEST_MULTIPLY) || defined(TEST_EXPONENTIATION) || defined(TEST_MULTIPLICATION) || defined(TEST_SIN) || defined(TEST_COS) || defined(TEST_EVALUATE) || defined(TEST_EXPONENTIATE) || defined(TEST_GET_BASE) || defined(TEST_GET_EXPONENT) || defined(TEST_IS_EXPONENTIATION) || defined(TEST_IS_MULTIPLICATION)
  using Exponent1 = Constant<3, 2>;
  using Symbol1 = Symbol<1>;
  using Power1 = Power<Symbol1, Exponent1>;
#endif

#ifdef TEST_POWER
  Power1 power;
#endif

#if defined(TEST_PRODUCT) || defined(TEST_SUM) || defined(TEST_MULTIPLY) || defined(TEST_EXPONENTIATION) || defined(TEST_MULTIPLICATION) || defined(TEST_SIN) || defined(TEST_COS) || defined(TEST_EVALUATE) || defined(TEST_EXPONENTIATE) || defined(TEST_GET_BASE) || defined(TEST_GET_EXPONENT) || defined(TEST_IS_EXPONENTIATION) || defined(TEST_IS_MULTIPLICATION)
  using Constant1 = Constant<42, 13>;
  using Symbol2 = Symbol<2>;
  using Symbol3 = Symbol<3>;
  using Power2 = Power<Symbol2, Exponent1>;
  using Power3 = Power<Symbol3, Exponent1>;
  using Product1 = Product<Constant1, Symbol1>;
  using Product2 = Product<Constant1, Power1>;
  using Product4 = Product<Symbol1, Symbol2>;
  using Product4a = Product<Symbol2, Symbol3>;
  using Product3 = Product<Constant1, Product4>;
  using Product5 = Product<Symbol1, Power2>;
  using Product6 = Product<Symbol1, Product4a>;
  using Product7 = Product<Power1, Symbol2>;
  using Product8 = Product<Power1, Power2>;
  using Product9 = Product<Power1, Product4a>;
#endif

#ifdef TEST_PRODUCT
  Product1 product1;
  Product2 product2;
  Product3 product3;
  Product4 product4;
  Product5 product5;
  Product6 product6;
  Product7 product7;
  Product8 product8;
  Product9 product9;
#endif

#if defined(TEST_SUM) || defined(TEST_EXPONENTIATION) || defined(TEST_MULTIPLICATION) || defined(TEST_SIN) || defined(TEST_COS) || defined(TEST_EVALUATE) || defined(TEST_EXPONENTIATE) || defined(TEST_GET_BASE) || defined(TEST_GET_EXPONENT) || defined(TEST_IS_EXPONENTIATION) || defined(TEST_IS_MULTIPLICATION)
  using Sum1 = Sum<Constant1, Symbol1>;
  using Sum2 = Sum<Constant1, Power1>;
  using Sum3 = Sum<Constant1, Product1>;
  using Sum4 = Sum<Constant1, Sum1>;
  using Sum5 = Sum<Symbol1, Symbol2>;
  using Sum6 = Sum<Symbol1, Power1>;
  using Sum7 = Sum<Symbol1, Product1>;
  using Sum8 = Sum<Symbol1, Sum1>;
  using Sum9 = Sum<Power1, Symbol2>;
  using Sum10 = Sum<Power1, Power2>;
  using Sum11 = Sum<Power1, Product1>;
  using Sum12 = Sum<Power1, Sum1>;
  using Sum13 = Sum<Product1, Symbol1>;
  using Sum14 = Sum<Product1, Power1>;
  using Sum15 = Sum<Product1, Product1>;
  using Sum16 = Sum<Product1, Sum1>;
#endif

#ifdef TEST_MULTIPLY
  multiply_t<Product8, Product9> foo;
#endif

#ifdef TEST_ADD_EQUALS
  add_equals<Symbol<0>, Symbol<0>> foo;
#endif

#ifdef TEST_SUM
  Sum1 sum1;
  Sum2 sum2;
  Sum3 sum3;
  Sum4 sum4;
  Sum5 sum5;
  Sum6 sum6;
  Sum7 sum7;
  Sum8 sum8;
  Sum9 sum9;
  Sum10 sum10;
  Sum11 sum11;
  Sum12 sum12;
  Sum13 sum13;
  Sum14 sum14;
  Sum15 sum15;
  Sum16 sum16;
#endif

#ifdef TEST_EXPONENTIATION
  Exponentiation<Sum1, Exponent1> exponentiation;
#endif

#ifdef TEST_MULTIPLICATION
  Multiplication<Sum1, Sum16> multiplication;
#endif

#if defined(TEST_SIN) || defined(TEST_COS)
  Sin<Symbol<1>> sin;
  Cos<Exponentiation<Sum1, Exponent1>> cos;
#endif

#ifdef TEST_INVERT
  invert_t<Symbol<1>> invert;
#endif

#ifdef TEST_NEGATE
  negate_t<Symbol<1>> negate;
#endif

#ifdef TEST_ADD
  add_t<Symbol<1>, Symbol<2>> add;
#endif

#ifdef TEST_EVALUATE
  Sum1 sum1;
  double val = evaluate(sum1);
#endif

#ifdef TEST_EXPONENTIATE
  exponentiate_t<Sum15, Exponent1> exponentiate;
#endif

#ifdef TEST_GET_BASE
  get_base_t<exponentiate_t<Sum15, Exponent1>> foo;
#endif

#ifdef TEST_GET_EXPONENT
  get_exponent_t<exponentiate_t<Sum15, Exponent1>> foo;
#endif

#ifdef TEST_GET_CONSTANT_FACTOR
  get_constant_factor_t<Product<Constant<42, 1>, Symbol<2>>> foo;
#endif

#ifdef TEST_GET_NON_CONSTANT_FACTOR
  get_nonconstant_factor_t<Product<Constant<42, 1>, Power<Symbol<2>, Constant<3, 2>>>> foo;
#endif

#ifdef TEST_IS_SYMBOL
  static_assert(is_symbol_v<Symbol<0>>, "fail");
#endif

#ifdef TEST_IS_POWER
  static_assert(is_power_v<Power<Symbol<0>, Constant<13, 42>>>, "fail");
#endif

#ifdef TEST_IS_PRODUCT
  static_assert(is_product_v<Product<Symbol<1>, Power<Symbol<1>, Constant<13, 42>>>>, "fail");
#endif

#ifdef TEST_IS_EXPONENTIATION
  static_assert(is_exponentiation_v<Exponentiation<Sum13, Exponent1>>, "fail");
#endif

#ifdef TEST_IS_MULTIPLICATION
  static_assert(is_multiplication_v<Multiplication<Sum1, Exponentiation<Sum13, Exponent1>>>, "fail");
#endif

#ifdef TEST_IS_SUM
  static_assert(is_sum_v<Sum<Symbol<0>, Symbol<1>>>, "fail");
#endif

#ifdef TEST_IS_SIN
  static_assert(is_sin_v<Sin<Symbol<1>>>, "fail");
#endif

#ifdef TEST_IS_COS
  static_assert(is_cos_v<Cos<Symbol<1>>>, "fail");
#endif
}
