#pragma once

#include <cstdint>
#include <array>

namespace symbolic {

// The following table shows when an operation needs parenthesis.
//
//            If we have a:
//              negation   |    power    |   product   |    ratio    |    sum      | difference
// of a:     +-----------------------------------------+----------------------------------------
// negation  |   -(-s)     |  (-s)^n     | x * -y      |   p/-q      | a + -b      | i - -j
// power     |   -(s^n)    |  (s^n)^m    | x * s^n     |   p/s^n     | a + s^n     | i - j^n       of_a_expr (vertical) is used
// product   |   -(x * y)  |  (x * y)^n  | x * y * z   |   x * y / q | a + x * y   | i - x * y   < before
//                                                         p/(x * y)                             < after
// ratio     |   -(p/q)    |  (p/q)^n    | x * p/q     |   p / q / o | a + p/q     | i - p/q     < before
//                                                         o/(p/q)                               < after
// sum       |   -(a + b)  |  (a + b)^n  | x * (a + b) |   p/(a + b) | a + b + c   | a + b - j   < before
//                                                                                   i - (a + b) < after
// difference|   -(i - j)  |  (i - j)^n  | x * (i - j) |   p/(i - j) | (i - j) + a | i - j - h   < before
//                                                     |             | a + i - j   | h - (i - j) < after
//                                                                                                 (binary operator) a_expr.

struct ParenthesisTableBase
{
  using mask_type = uint32_t;   // [difference][sum][ratio][product][power][negation][symbol]  <-- horizontal
                                // where each [operation] is: {after}{before}
                                // where each {position} is: 0 or 1 (one bit).
  enum Operation { symbol, negation, power, product, ratio, sum, difference, operation_count };
  static constexpr int position_width = 2;                                      // Width of {position} in bits.
  static constexpr int operation_width = position_width * operation_count;      // Width of [operation] in bits.

  static consteval mask_type before(Operation op) { return mask_type{1} << (op * position_width); }
  static consteval mask_type after(Operation op) { return mask_type{2} << (op * position_width); }
  static consteval mask_type with(Operation op) { return mask_type{3} << (op * position_width); }      // Both before and after.
};

struct ParenthesisTable : ParenthesisTableBase
{
  static constexpr std::array<uint32_t, operation_count> table = {
    //           when used with: negation         power         product         ratio          sum           difference
    // This needs -.
    // parenthesis v
             /*symbol*/     0, // Never
             /*negation*/   with(negation) | with(power),
             /*power*/      with(negation) | with(power),
             /*product*/    with(negation) | with(power) |                after(ratio),
             /*ratio*/      with(negation) | with(power) |                after(ratio),
             /*sum*/        with(negation) | with(power) | with(product) | with(ratio) |               after(difference),
             /*difference*/ with(negation) | with(power) | with(product) | with(ratio) | before(sum) | after(difference)
  };
  // encoded returns a double mask: <vertical><horizontal>, where <horizontal> is the mask showns above
  // and <vertical> just encodes the op as a mask: have both bits corresponding to the op set.
  static consteval mask_type encoded(Operation op) { return (with(op) << operation_width) | table[op]; }
};

enum class precedence : ParenthesisTable::mask_type
{
  constant  = 0,
  symbol     = constant,
  negation   = ParenthesisTable::encoded(ParenthesisTable::negation),
  power      = ParenthesisTable::encoded(ParenthesisTable::power),
  product    = ParenthesisTable::encoded(ParenthesisTable::product),
  ratio      = ParenthesisTable::encoded(ParenthesisTable::ratio),
  sum        = ParenthesisTable::encoded(ParenthesisTable::sum),
  difference = ParenthesisTable::encoded(ParenthesisTable::difference)
};

enum before_or_after : ParenthesisTable::mask_type
{
  before = 0x1555,
  after = 0x2aaa
};

// Return true when `of_a_expr` needs parenthesis when used before/after `a_expr`
// (in case `a_expr` is a binary operator; if it is not then position should be `before`).
constexpr inline bool needs_parens(precedence of_a_expr, precedence a_expr, before_or_after position)
{
  return (static_cast<ParenthesisTable::mask_type>(of_a_expr) >> ParenthesisTable::operation_width) &
    static_cast<ParenthesisTable::mask_type>(a_expr) & position;
}

} // namespace symbolic
