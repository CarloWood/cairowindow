#pragma once

#include "IdRange.h"
#include "Power.h"
#include "Expression.h"

namespace symbolic {

template<Expression E1, Expression E2>
class Product : public ExpressionTag
{
 public:
  static constexpr precedence s_precedence = precedence::product;
  static constexpr IdRange<E1::id_range.begin, E2::id_range.end> id_range;

 private:
  E1 arg1_;
  E2 arg2_;

 public:
  inline constexpr Product(E1 const& arg1, E2 const& arg2);

  constexpr E1 const& arg1() const { return arg1_; }
  constexpr E2 const& arg2() const { return arg2_; }

#ifdef SYMBOLIC_PRINTING
  constexpr bool needs_parens(before_or_after position, precedence prec) const { return prec <= s_precedence; }
  constexpr bool needs_parens(precedence prec) const { return prec <= s_precedence; }

  void print_on(std::ostream& os) const
  {
    bool need_parens = arg1_.needs_parens(before, s_precedence);
    if (need_parens)
      os << '(';
    arg1_.print_on(os);
    if (need_parens)
      os << ')';
    os << " * ";
    need_parens = arg2_.needs_parens(after, s_precedence);
    if (need_parens)
      os << '(';
    arg2_.print_on(os);
    if (need_parens)
      os << ')';
  }
#endif
};

template<typename T>
struct is_product : std::false_type { };

template<typename T>
constexpr bool is_product_v = is_product<T>::value;

template<Expression E1, Expression E2>
struct is_product<Product<E1, E2>> : std::true_type { };

template<typename E>
concept ProductType = is_product_v<E>;

template<Expression E1, Expression E2>
constexpr Product<E1, E2>::Product(E1 const& arg1, E2 const& arg2) : arg1_(arg1), arg2_(arg2)
{
  static_assert(!is_product_v<E1>, "The first factor of a Product must not be a Product itself.");
}

template<Expression E1, Expression E2, Expression E3, Expression E4>
constexpr bool is_same_expression(Product<E1, E2> const& arg1, Product<E3, E4> const& arg2)
{
  return std::is_same_v<E1, E3> && std::is_same_v<E2, E4>;
}

// The design of the following multiplication operator is as follows:
//
// Prerequisite is that all Product have the form {non-Product, Expression}
// and that they are already ordered; in other words, each Product has
// a cannonical form: if they represent the same product, their types
// will be the same.
//
// This also means that every Product can be thought of as a chaining
// of non-Product's, because if in the above Expression is a Product,
// it itself begins again with a non-Product.
//
// Lets call the non-Products involved, a, b, c, d etc. Then a Product
// can be thought of as an ordered series: b * (c * (f * (h * x))).
//
// For Symbols the ordering is determined by their Id. There can only
// be one Constant because those contract eachother; a Constant is
// always "less" than a Symbol with respect to product ordering.
// Hence, in the above (only) 'b' can be a Constant, but that is not
// necessary.
//
// The ordering is determined by the smallest Symbol Id in the Expression
// (or -1 for a Constant). However, each Expression also references the
// largest Symbol/Constant Id (plus one) as `end`. This allows us to quickly
// determine if a product has an overlapping range or not.
//
// In the comments below I use '*' for an operator* operation and an 'x'
// for a Product operation. Thus (a x b x c) really is (a x (b x c))
// (or rather, Product{a, Product{b, c}}).
//
template<Expression E1, Expression E2>
constexpr auto operator*(E1 const& arg1, E2 const& arg2)
{
  if constexpr (E1::id_range < E2::id_range)                    // e.g. (a x b x c) * (d x e x f), or a * (d x e x f)   [non-overlapping ranges].
  {
    if constexpr (is_product_v<E1>)
      return Product{arg1.arg1(), arg1.arg2() * arg2};          // a x ((b x c) * (d x e x f)).
    else
      return Product{arg1, arg2};                               // arg1 is a non-Product; keep the order (aka a x (d x e x f)).
  }
  else if constexpr (E2::id_range < E1::id_range)               // e.g. (d x e x f) * (a x b x c), or (d x e x f) * a  [non-overlapping ranges].
  {
    if constexpr (is_product_v<E2>)
      return Product{arg2.arg1(), arg2.arg2() * arg1};          // a x ((b x c) * (d x e x f)).
    else
      return Product{arg2, arg1};                               // arg2 is a non-Product; swap the order (aka a x (d x e x f)).
  }
  // The ranges overlap.
  else if constexpr (E2::id_range.begin < E1::id_range.begin)   // Make sure that E1::id_range.begin is not larger than that of E2.
    return arg2 * arg1;
  else if constexpr (is_product_v<E1> && is_symbol_v<E2>)       // e.g. (a x b x d) * c, where c >= a.
    return arg1.arg1() * (arg2 * arg1.arg2());                  // a * (c * (b x d)).
  else
  {
    static_assert(!(is_product_v<E1> && is_product_v<E2>), ""); // For this case the specialization below is called.
    ASSERT(is_same_expression(arg1, arg2));
    return Power<E1, 2, 1>{arg1};
  }
}

// This is a work-around for the fact that you can't recursively call
// a template function with a auto-deduced return type "before it is
// defined".
template<Expression E1, Expression E2, Expression E3, Expression E4>
requires (!(Product<E1, E2>::id_range < Product<E3, E4>::id_range) &&           // The range of the two products
          !(Product<E3, E4>::id_range < Product<E1, E2>::id_range) &&           // must overlap.
          !(Product<E3, E4>::id_range.begin < Product<E1, E2>::id_range.begin)) // The id of the start of the range of the first product
                                                                                // maybe not be larger than the id of the start of the
                                                                                // range of the second product.
constexpr auto operator*(Product<E1, E2> const& arg1, Product<E3, E4> const& arg2)
{
  // e.g. (a x c x e) * (b x d x f), where b >= a --> a has the smallest id.
  if constexpr (is_product_v<E2>)
    return Product{arg1.arg1(), arg2.arg1() * (arg1.arg2() * arg2.arg2())};
  else
  return Product{arg1.arg1(), arg2 * arg1.arg2()};              // a x ((b x d x f) * (c x e))
}

} // namespace symbolic
