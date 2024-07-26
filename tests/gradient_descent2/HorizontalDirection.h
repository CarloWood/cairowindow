#pragma once

#include <string>
#include "debug.h"
#ifdef CWDEBUG
#include <ostream>
#endif

namespace gradient_descent {

// Are we looking for the next extreme left or right of the previous extreme found?
enum class HorizontalDirection
{
  left = -1,
  undecided = 0,
  right = 1,
};

// The output of find_extreme:
//   - the region where the requested (or returned) extreme is found in.
//
//      ------------|---------------|--------------
//                  wl              wr
// output:          |               |
//      ---left---->|               |
//                  |<--inbetween-->|
//                  |               |<----right----
enum class Region
{
  left = -1,
  inbetween = 0,
  right = 1,
  unknown = 2,
  invalid = 3
};

inline HorizontalDirection opposite(HorizontalDirection hdirection)
{
  return static_cast<HorizontalDirection>(-static_cast<int>(hdirection));
}

// Return if x is left or right of y.
inline HorizontalDirection side(double x, double y)
{
  ASSERT(x != y);
  return x < y ? HorizontalDirection::left : HorizontalDirection::right;
}

// Return if step is going left or right.
inline HorizontalDirection side(double step)
{
  ASSERT(step != 0.0);
  return step < 0.0 ? HorizontalDirection::left : HorizontalDirection::right;
}

// Add -/+ delta to x so that we go left/right as per hdirection.
inline double add_to(double x, double delta, HorizontalDirection hdirection)
{
  ASSERT(delta > 0.0);
  return x + static_cast<int>(hdirection) * delta;
}

class HorizontalDirectionToInt {
 private:
  int val_;
 public:
  HorizontalDirectionToInt(HorizontalDirection hdirection) : val_(static_cast<int>(hdirection)) { }
  HorizontalDirectionToInt(Region region) : val_(static_cast<int>(region)) { ASSERT(val_ == -1 || val_ == 1); }
  operator int() const { return val_; }
  int as_index() const { return (val_ + 1) >> 1; }
};

std::string to_string(HorizontalDirection hdirection);
std::string to_string(Region region);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, HorizontalDirection hdirection)
{
  return os << to_string(hdirection);
}

inline std::ostream& operator<<(std::ostream& os, Region region)
{
  return os << to_string(region);
}
#endif

} // namespace gradient_descent
