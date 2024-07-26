#pragma once

#include <string>
#ifdef CWDEBUG
#include <iostream>
#endif
#include "debug.h"

namespace gradient_descent {

enum class CubicToNextSampleType : int
{
  unknown = 0,  // ?? Used when SampleNode::next_ is nullptr.
  // These types have no minimum or maximum:
  up,           // /  This and the next sample have a positive derivative and the cubic has no extremes in between.
  down,         // \  This and the next sample have a negative derivative and the cubic has no extremes in between.
  right_stop,   // /^ This has a positive derivative and no extremes could be found at all on the right (up till the max. energy).
  left_stop,    // ^\ The next sample has a negative derivative and no extremes could be found at all on the left (up till the max. energy).
  right_min,    // _/ This sample is a local minimum; the next sample has a positive derivative.
  left_min,     // \_ The next sample is a local minimum; this sample has a negative derivative.
  right_max,    // ‾\ This sample is a local maximum; the next sample has a negative derivative.
  left_max,     // /‾ The next sample is a local maxmimum; this sample has a positive derivative.
  // These types have a minimum:
  min,          // \/ This sample has negative derivative and the next sample has a positive derivative.
                //    The cubic has a minimum in between, but no maximum.
  // These types have a maximum (too):
  min_max,      // \/\ This and the next sample have a negative derivative and the cubic has both derivative in between.
  max_min,      // /\/ This and the next sample have a positive derivative and the cubic has both derivative in between.
  // This type doesn't have a minimum, only a maximum:
  max,          // /\ This sample has a positive derivative and the next sample has a negative derivative.
                // The cubic has a maximum in between, but no minimum.
};

inline bool has_minimum(CubicToNextSampleType type)
{
  // We can't answer this question if type is unknown.
  ASSERT(type != CubicToNextSampleType::unknown);
  return CubicToNextSampleType::min <= type && type < CubicToNextSampleType::max;
}

inline bool has_maximum(CubicToNextSampleType type)
{
  // We can't answer this question if type is unknown.
  ASSERT(type != CubicToNextSampleType::unknown);
  return type > CubicToNextSampleType::min;
}

std::string to_string(CubicToNextSampleType state);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, CubicToNextSampleType state)
{
  return os << to_string(state);
}
#endif

} // namespace gradient_descent
