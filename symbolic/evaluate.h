#pragma once

#include "expression_traits.h"

namespace symbolic {

template<Expression T>
inline double evaluate(T const&);

template<ConstantType T>
double evaluate()
{
  return static_cast<double>(T::s_enumerator) / T::s_denominator;
}

template<SymbolType T>
double evaluate()
{
  return T::get_value(T::s_id);
}

template<PowerType T>
double evaluate()
{
  return std::pow(evaluate<typename T::base_type>(), evaluate<typename T::exponent_type>());
}

template<ProductType T>
double evaluate()
{
  return evaluate<typename T::arg1_type>() * evaluate<typename T::arg2_type>();
}

template<ExponentiationType T>
double evaluate();

template<MultiplicationType T>
double evaluate();

template<SumType T>
double evaluate()
{
  return evaluate<typename T::arg1_type>() + evaluate<typename T::arg2_type>();
}

template<ExponentiationType T>
double evaluate()
{
  return std::pow(evaluate<typename T::base_type>(), evaluate<typename T::exponent_type>());
}

template<MultiplicationType T>
double evaluate()
{
  return evaluate<typename T::arg1_type>() * evaluate<typename T::arg2_type>();
}

template<Expression T>
double evaluate(T const&)
{
  return evaluate<T>();
}

} // namespace symbolic
