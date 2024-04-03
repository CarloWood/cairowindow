#pragma once

#include "Expression.h"
#include "is_constant.h"
#include "is_power.h"
#include "is_product.h"
#include "get_base.h"

namespace symbolic {

template<Expression E1, Expression E2>
struct is_less_Product : std::false_type { };

template<Expression E1, Expression E2>
constexpr bool is_less_Product_v = is_less_Product<E1, E2>::value;

template<int Enumerator1, int Denominator1, Expression E2>
requires (!is_constant_v<E2>)
struct is_less_Product<Constant<Enumerator1, Denominator1>, E2> : std::true_type { };

template<int Id1, int Id2>
requires (Id1 < Id2)
struct is_less_Product<Symbol<Id1>, Symbol<Id2>> : std::true_type { };

template<int Id1, int Id2, ConstantType Exponent>
requires (Id1 < Id2)
struct is_less_Product<Symbol<Id1>, Power<Symbol<Id2>, Exponent>> : std::true_type { };

template<int Id1, ConstantType Exponent1, int Id2>
requires (Id1 < Id2)
struct is_less_Product<Power<Symbol<Id1>, Exponent1>, Symbol<Id2>> : std::true_type { };

template<int Id1, ConstantType Exponent1, int Id2, ConstantType Exponent2>
requires (Id1 < Id2)
struct is_less_Product<Power<Symbol<Id1>, Exponent1>, Power<Symbol<Id2>, Exponent2>> : std::true_type { };

template<int Id, Expression E1, Expression E2>
requires (Id < get_base_t<E1>::s_id)
struct is_less_Product<Symbol<Id>, Product<E1, E2>> : std::true_type { };

template<int Id, ConstantType Exponent, Expression E1, Expression E2>
requires (Id < get_base_t<E1>::s_id)
struct is_less_Product<Power<Symbol<Id>, Exponent>, Product<E1, E2>> : std::true_type { };

template<Expression E1, Expression E2, int Id>
requires (get_base_t<E1>::s_id < Id)
struct is_less_Product<Product<E1, E2>, Symbol<Id>> : std::true_type { };

template<Expression E1, Expression E2, int Id, ConstantType Exponent>
requires (get_base_t<E1>::s_id < Id)
struct is_less_Product<Product<E1, E2>, Power<Symbol<Id>, Exponent>> : std::true_type { };

template<Expression E1, Expression E2, Expression E3, Expression E4>
requires (get_base_t<E1>::s_id < get_base_t<E3>::s_id ||
    (get_base_t<E1>::s_id == get_base_t<E3>::s_id && get_base_t<E2>::s_id < get_base_t<E4>::s_id))
struct is_less_Product<Product<E1, E2>, Product<E3, E4>> : std::true_type { };

} // namespace symbolic
