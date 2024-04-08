#pragma once

#include "expression_traits.h"
#include <cmath>

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

// Forward declarations of all evaluate specializations.

template<ExponentiationType T>
double evaluate();

template<MultiplicationType T>
double evaluate();

template<SinType T>
double evaluate();

template<CosType T>
double evaluate();

template<LogType T>
double evaluate();

template<AtanType T>
double evaluate();

//

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

template<SinType T>
double evaluate()
{
  return std::sin(evaluate<typename T::arg_type>());
}

template<CosType T>
double evaluate()
{
  return std::cos(evaluate<typename T::arg_type>());
}

template<LogType T>
double evaluate()
{
  return std::log(evaluate<typename T::arg_type>());
}

template<AtanType T>
double evaluate()
{
  return std::atan(evaluate<typename T::arg_type>());
}

} // namespace symbolic
