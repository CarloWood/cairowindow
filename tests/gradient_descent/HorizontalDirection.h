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

// In the context of find_extreme, two samples are involved -
// the left-most sample (at wl) and the right-most sample (at wr).
//
// As input of find_extreme:
//   - an optional restriction on the region(s) where the returned extreme is allowed to be in.
//
//      ------------|---------------|--------------
//                  wl              wr
// input:           |               |
//      -----------left------------>|                   : only returns inbetween or left.
//                  |<--------------right----------     : only returns inbetween or right.
//      -------------------none--------------------     : returns whatever.
enum class Restriction
{
  left = -1,
  none = 0,
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
  unknown = 2
};

inline bool operator==(Region region, Restriction restriction)
{
  // Restriction; Region: -1   0     1
  //     -1             | T  | T  | F  |
  //      0             | T  | T  | T  |
  //      1             | F  | T  | T  |
  return static_cast<int>(region) * static_cast<int>(restriction) == -1;
}

inline HorizontalDirection opposite(HorizontalDirection hdirection)
{
  return static_cast<HorizontalDirection>(-static_cast<int>(hdirection));
}

class HorizontalDirectionToInt {
 private:
  int val_;
 public:
  HorizontalDirectionToInt(HorizontalDirection hdirection) : val_(static_cast<int>(hdirection)) { ASSERT(val_ == -1 || val_ == 1); }
  HorizontalDirectionToInt(Restriction restriction) : val_(static_cast<int>(restriction)) { ASSERT(val_ == -1 || val_ == 1); }
  HorizontalDirectionToInt(Region region) : val_(static_cast<int>(region)) { ASSERT(val_ == -1 || val_ == 1); }
  operator int() const { return val_; }
};

std::string to_string(HorizontalDirection hdirection);
std::string to_string(Region region);
std::string to_string(Restriction region);

#ifdef CWDEBUG
inline std::ostream& operator<<(std::ostream& os, HorizontalDirection hdirection)
{
  return os << to_string(hdirection);
}

inline std::ostream& operator<<(std::ostream& os, Region region)
{
  return os << to_string(region);
}

inline std::ostream& operator<<(std::ostream& os, Restriction restriction)
{
  return os << to_string(restriction);
}
#endif

} // namespace gradient_descent
