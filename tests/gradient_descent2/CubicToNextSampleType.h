#pragma once

#include <string>
#include "ExtremeType.h"
#include "HorizontalDirection.h"
#ifdef CWDEBUG
#include <iostream>
#endif
#include "debug.h"

namespace gradient_descent {

static constexpr int minimum_bit   = 0x01;
static constexpr int maximum_bit   = 0x02;
static constexpr int rising_bit    = 0x04;
static constexpr int falling_bit   = 0x08;
static constexpr int left_min_bit  = 0x10;
static constexpr int left_max_bit  = 0x20;
static constexpr int right_min_bit = 0x40;
static constexpr int right_max_bit = 0x80;
static constexpr int max_before_min_bit = 0x100;
static constexpr int id_shift = 9;

enum class CubicToNextSampleType : int
{
  unknown = 0,  // ?? Used when SampleNode::next_ is nullptr.

  // These types have no minimum or maximum:

  // __  This and the next sample have a derivative of zero and the cubic is a straight horizontal line.
  flat                  =  1 << id_shift,

  // /   This and the next sample have a positive derivative and the cubic has no extremes in between.
  up                    =  2 << id_shift | rising_bit,

  // \   This and the next sample have a negative derivative and the cubic has no extremes in between.
  down                  =  3 << id_shift | falling_bit,

  // /^  This has a positive derivative and no extremes could be found at all on the right (up till the max. energy).
  right_stop            =  4 << id_shift | rising_bit,

  // ^\  The next sample has a negative derivative and no extremes could be found at all on the left (up till the max. energy).
  left_stop             =  5 << id_shift | falling_bit,

  // _/  This sample is a local minimum; the next sample has a positive derivative.
  right_min             =  6 << id_shift | rising_bit | right_min_bit,

  // \_  The next sample is a local minimum; this sample has a negative derivative.
  left_min              =  7 << id_shift | falling_bit | left_min_bit,

  // ‾\  This sample is a local maximum; the next sample has a negative derivative.
  right_max             =  8 << id_shift | falling_bit | right_max_bit,

  // /‾  The next sample is a local maxmimum; this sample has a positive derivative.
  left_max              =  9 << id_shift | rising_bit | left_max_bit,

  // ‾\_ This sample is a local maximum; the next sample is a local minimum.
  right_max_left_min    = 10 << id_shift | falling_bit | right_max_bit | left_min_bit | max_before_min_bit,

  // _/‾ This sample is a local minimum; the next sample is a local maximum.
  right_min_left_max    = 11 << id_shift | rising_bit | right_min_bit | left_max_bit,

  // These types have a minimum:

  // \/  This sample has negative derivative and the next sample has a positive derivative.
  min                   = 12 << id_shift | minimum_bit,

  // ‾\/ This sample is a local maximum and the next sample has a positive derivative.
  right_max_min         = 13 << id_shift | minimum_bit | right_max_bit | max_before_min_bit,

  // \/‾ The next sample is a local maximum and this sample has a negative derivative.
  min_left_max          = 14 << id_shift | minimum_bit | left_max_bit,

  // These types have a maximum (too):

  // \/\ This and the next sample have a negative derivative and the cubic has both derivative in between.
  min_max               = 15 << id_shift | minimum_bit | maximum_bit,

  // /\/ This and the next sample have a positive derivative and the cubic has both derivative in between.
  max_min               = 16 << id_shift | minimum_bit | maximum_bit | max_before_min_bit,

  // This type doesn't have a minimum, only a maximum:

  // /\  This sample has a positive derivative and the next sample has a negative derivative.
  max                   = 17 << id_shift | maximum_bit,

  // /\_ The next sample is a local minimum and this sample has a positive derivative.
  max_left_min          = 18 << id_shift | maximum_bit | left_min_bit | max_before_min_bit,

  // _/\ This sample is a local minimum and the next sample has a negative derivative.
  right_min_max         = 19 << id_shift | maximum_bit | right_min_bit,
};

inline bool has_unfound_minimum(CubicToNextSampleType type)
{
  // We can't answer this question if type is unknown.
  ASSERT(type != CubicToNextSampleType::unknown);
  return static_cast<int>(type) & minimum_bit;
}

inline bool has_unfound_maximum(CubicToNextSampleType type)
{
  // We can't answer this question if type is unknown.
  ASSERT(type != CubicToNextSampleType::unknown);
  return static_cast<int>(type) & maximum_bit;
}

inline bool has_minimum(CubicToNextSampleType type)
{
  // We can't answer this question if type is unknown.
  ASSERT(type != CubicToNextSampleType::unknown);
  return static_cast<int>(type) & (minimum_bit|left_min_bit|right_min_bit);
}

inline bool has_maximum(CubicToNextSampleType type)
{
  // We can't answer this question if type is unknown.
  ASSERT(type != CubicToNextSampleType::unknown);
  return static_cast<int>(type) & (maximum_bit|left_max_bit|right_max_bit);
}

// extreme_type <---direction--- opposite extreme_type.
inline bool has_extreme_order(CubicToNextSampleType type, ExtremeType extreme_type, HorizontalDirection direction)
{
  ASSERT(extreme_type != ExtremeType::unknown);

  int b = (static_cast<int>(type) & max_before_min_bit) ? 0 : 2;        // 0: max_min, 2: min_max.
  int e = static_cast<int>(extreme_type) + 1;                           // 0: minimum, 2: maximum.
  int d = static_cast<int>(direction) + 1;                              // 0: left,    2: right.
  bool correct_order = b ^ d ^ e;

  return has_maximum(type) && has_minimum(type) && correct_order;
}

std::string to_string(CubicToNextSampleType type);
char const* to_utf8_art(CubicToNextSampleType type);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, CubicToNextSampleType type)
{
  return os << to_string(type);
}
#endif

} // namespace gradient_descent
