#include "sys.h"
#include "Scale.h"
#include "utils/macros.h"

namespace gradient_descent {

double Scale::calculate_value(double const inflection_point_w) const
{
  // The (value of the) scale is losely defined at the distance over which w can change
  // without L(w) starting to deviate significantly from the cubic.
  //
  // Moreover, because of the nature of the algorithm, which approaches a critical point
  // from one side, we don't want the scale to be much different if one side hasn't be
  // explored yet.
  //
  // For example,
  //
  //      cp          : critical point
  //     l|           : left-most sample
  //     ↓↓
  //     .-. r        : right-most sample
  //    /   \↓
  //   /     \
  //
  // If l is very close to cp, or non-existant, which gives a distance of 0, then
  // we just want scale to be `r - cp`; assuming that there is no reason that
  // L(w) would look different on the left side of cp if we never explored that side.
  //
  // But if `cp - l` is significant, although less than `r - cp` the chance increases
  // that we DID explore the left side and there is too much deviation left of l.
  //
  // Therefore, we return a "weighted" average between `cp - l` and `r - cp`, weighted
  // with their own value: (w1*A + w2*B) / (w1 + w2), where we'll use w1=A and w2=B.
  // In this case A = cp - l and B = r - cp. Hence,
  //
  //   value = (A^2 + B^2) / (A + B)
  //
  // Where A and B are positive if l and r are on opposite sides of cp.
  //
  // However, if l and r are on the same side of cp, for example,
  //
  //      cp          : critical point
  //     r|           : right-most sample (now B is negative)
  //     ↓↓
  //   l .-.          : left-most sample
  //   ↓/   \
  //   /     \
  //
  // Then we want the returned scale to be `r - l`.
  //
  // Extension with inflection point
  // -------------------------------
  //
  // Let C be the critical point at critical_point_w_.
  // Let I be the inflection point of the cubic.
  // Then the following situations are possible:
  //
  //                 C/I                                    : C and I on top of eachother (C is the inflection point)
  //                                                          this happens when the cubic has no local extrema.
  //
  //           I     C                                      : The inflection point is on the left of the critical point.
  //                                                          \      C                    /
  //                                                           \    /\             /\    /
  //                                                            \  I  \     or    /  I  /
  //                                                             \/    \         /    \/
  //                                                                    \       /      C
  //                 C     I                                : The inflection point is on the right of the critical point.
  //                                                          \                    C      /
  //                                                           \    /\             /\    /
  //                                                            \  I  \     or    /  I  /
  //                                                             \/    \         /    \/
  //                                                             C      \       /
  //
  // Let the left_edge_w_ be l and the right_edge_w_ sample be r.
  // Placing l and r in each of the above cases we have the possibilities:
  //
  //                                      Class (value returned by classify(l, I, C, r)).
  //   l < r < C = I   : r - l                4
  //   C = I < l < r   : r - l                4
  //   l < r < I < C   : r - l                4
  //   I < l < r < C   : r - l                4
  //   I < C < l < r   : r - l                4
  //   l < r < C < I   : r - l                4
  //   C < l < r < I   : r - l                4
  //   C < I < l < r   : r - l                4
  //   l < I < r < C   : wa(l, I, r)          1
  //   C < l < I < r   : wa(l, I, r)          1
  //   l < C = I < r   : wa(l, C, r)          3 - however, classify returns 2.
  //   I < l < C < r   : wa(l, C, r)          3
  //   l < C < r < I   : wa(l, C, r)          3
  //   l < I < C < r   : wa(I, C, r)          0
  //   l < C < I < r   : wa(l, C, I)          2
  //
  // where wa(l, C, r) is the weighted average around C: ((C - l)^2 + (r - C)^2) / (r - l) etc.
  //
  // Note that the Class values have been chosen such that we can simply replace element `Class`
  // in the LCRI array with the inflection point, replacing respectively l, C or r (or I with itself
  // in the case of Class 3) and then call weighted_average with that altered array.

  DoutEntering(dc::notice, "SampleNode::calculate_value(" << inflection_point_w << ")");
  Dout(dc::notice, "type = " << type_ << "; critical_point_w = " << critical_point_w_ <<
      "; left_edge_w = " << left_edge_w_ << "; right_edge_w = " << right_edge_w_);
  // critical_point_w_ should be initialized.
  ASSERT(type_ != CriticalPointType::none);
  LCRI_type LCRI = {left_edge_w_, critical_point_w_, right_edge_w_, inflection_point_w};
  int Class = classify(LCRI);
  Dout(dc::notice, "Class = " << Class);
  double result;
  if (Class == 4)
    result = LCRI[Ri] - LCRI[Li];
  else
  {
    if (AI_LIKELY(type_ != CriticalPointType::inflection_point))
      LCRI[Class] = LCRI[Ii];
    result = weighted_average(LCRI);
  }
  Dout(dc::notice, "Scale value: " << result);
  return result;
}

} // namespace gradient_descent
