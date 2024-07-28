#include "sys.h"
#include "SampleNode.h"
#include "AlgorithmEventType.h"

namespace gradient_descent {

void SampleNode::initialize_cubic(SampleNode const& next
    COMMA_CWDEBUG_ONLY(events::Server<AlgorithmEventType>& event_server, bool this_is_last)) const
{
  DoutEntering(dc::notice, "SampleNode::initialize_cubic(" << next << "..., " << std::boolalpha << this_is_last << ")");

  // next must always be on the right.
  ASSERT(w() < next.w());

  cubic_.initialize(w(), Lw(), dLdw(), next.w(), next.Lw(), next.dLdw());

#ifdef CWDEBUG
  // If the current object (this) was the last sample that was added, then the cubic is on the right
  // of the last sample. We use a different color for that. Otherwise assume next was the last.
  HorizontalDirection side = this_is_last ? HorizontalDirection::right : HorizontalDirection::left;
  event_server.trigger(AlgorithmEventType{cubic_polynomial_event, cubic_, side});
#endif

  // Get the sign of the derivatives.
  int const sign_dLdw_0 = dLdw() < 0.0 ? -1 : dLdw() > 0.0 ? 1 : 0;
  int const sign_dLdw_1 = next.dLdw() < 0.0 ? -1 : next.dLdw() > 0.0 ? 1 : 0;
  bool const neither_derivative_is_zero = (sign_dLdw_0 & sign_dLdw_1) != 0;

  using enum CubicToNextSampleType;

  // The easiest first.
  if (sign_dLdw_0 != sign_dLdw_1)
  {
    // This must be /\, \/ -- or, very unlikely, one of the derivatives is zero.
    if (AI_LIKELY(neither_derivative_is_zero))
    {
      type_ = (sign_dLdw_0 == 1)
          ? max                         // /\.
          : min;                        // \/
    }
    else
    {
      // One of the derivatives is zero!
      if (sign_dLdw_0 == 0)             // The right part of an extreme.
      {
        double d2Ldw2 = cubic_.second_derivative(w());
        if (sign_dLdw_1 == 1)
        {
          type_ = (d2Ldw2 > 0.0)
              ? right_min               // _/
              : right_max_min;          // ‾\/
        }
        else
        {
          type_ = (d2Ldw2 > 0.0)
              ? right_min_max           // _/\.
              : right_max;              // ‾\.
        }
      }
      else                              // The left part of an extreme.
      {
        double d2Ldw2 = cubic_.second_derivative(next.w());
        if (sign_dLdw_0 == 1)
        {
          type_ = (d2Ldw2 > 0.0)
              ? max_left_min          // /\_
              : left_max;             // /‾
        }
        else
        {
          type_ = (d2Ldw2 > 0.0)
              ? left_min              // \_
              : min_left_max;         // \/‾
        }
      }
    }
  }
  else
  {
    // The signs are equal, so this must be /, \, /\/, \/\ -- or, very unlikely, both derivatives are zero.
    if (AI_LIKELY(neither_derivative_is_zero))
    {
      // The signs are equal (and non-zero). Now we have to calculate the extremes of the cubic.
      // Put the minimum (or inflection point) in index 0.
      std::array<double, 2> extremes;
      int number_of_extremes = cubic_.get_extremes(extremes);

      if (number_of_extremes < 2 || extremes[0] < w() || next.w() < extremes[0])
      {
        // If only one extreme would fall in between the samples, then the derivatives can't have the same sign.
        ASSERT(number_of_extremes < 2 || extremes[1] < w() || next.w() < extremes[1]);
        // There are no extremes in between the samples.
        type_ = (sign_dLdw_0 == 1)
            ? up                        // /
            : down;                     // \.
      }
      else
      {
        // Both extremes are between the two samples.
        ASSERT(w() < extremes[1] && extremes[1] < next.w());
        type_ = (sign_dLdw_0 == 1)
            ? max_min                   // /\/
            : min_max;                  // \/\.
      }
    }
    else
    {
      // Both derivatives are zero!
      if (Lw() > next.Lw())             // ‾\_
        type_ = right_max_left_min;
      else if (Lw() < next.Lw())        // _/‾
        type_ = right_min_left_max;
      else
        type_ = flat;                   // __
    }
  }

  Dout(dc::notice, "Result: " << *this);
}

} // namespace gradient_descent
