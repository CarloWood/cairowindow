#pragma once

#include "canonical_constant.h"
#include "multiply_fwd.h"
#include "Constant.h"

namespace symbolic {

template<Expression E>
struct negate
{
  using type = typename multiply<Constant<-1, 1>, E>::type;
};

template<Expression E>
using negate_t = typename negate<E>::type;

template<int Enumerator, int Denominator>
struct negate<Constant<Enumerator, Denominator>>
{
  using type = canonical_constant_t<-Enumerator, Denominator>;
};

} // namespace symbolic
