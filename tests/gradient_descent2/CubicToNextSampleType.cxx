#include "sys.h"
#include "CubicToNextSampleType.h"
#include "utils/macros.h"
#include "debug.h"

namespace gradient_descent {

std::string to_string(CubicToNextSampleType type)
{
  using enum CubicToNextSampleType;
  switch (type)
  {
    AI_CASE_RETURN(unknown);
    AI_CASE_RETURN(flat);
    AI_CASE_RETURN(up);
    AI_CASE_RETURN(down);
    AI_CASE_RETURN(right_stop);
    AI_CASE_RETURN(left_stop);
    AI_CASE_RETURN(right_min);
    AI_CASE_RETURN(left_min);
    AI_CASE_RETURN(right_max);
    AI_CASE_RETURN(left_max);
    AI_CASE_RETURN(right_max_left_min);
    AI_CASE_RETURN(right_min_left_max);
    AI_CASE_RETURN(min);
    AI_CASE_RETURN(right_max_min);
    AI_CASE_RETURN(min_left_max);
    AI_CASE_RETURN(min_max);
    AI_CASE_RETURN(max_min);
    AI_CASE_RETURN(max);
    AI_CASE_RETURN(max_left_min);
    AI_CASE_RETURN(right_min_max);
  }
  AI_NEVER_REACHED
}

// For now assume the ostream will understand an UTF8 byte stream.
char const* to_utf8_art(CubicToNextSampleType type)
{
  using enum CubicToNextSampleType;
  switch (type)
  {
    case unknown:
      return "<unknown>";
    case flat:
      return "__";
    case up:
      return "/";
    case down:
      return "\\";
    case right_stop:
      return "/^";
    case left_stop:
      return "^\\";
    case right_min:
      return "_/";
    case left_min:
      return "\\_";
    case right_max:
      return reinterpret_cast<char const*>(u8"‾\\");
    case left_max:
      return reinterpret_cast<char const*>(u8"/‾");
    case right_max_left_min:
      return reinterpret_cast<char const*>(u8"‾\\_");
    case right_min_left_max:
      return reinterpret_cast<char const*>(u8"_/‾");
    case min:
      return "\\/";
    case right_max_min:
      return reinterpret_cast<char const*>(u8"‾\\/");
    case min_left_max:
      return reinterpret_cast<char const*>(u8"\\/‾");
    case min_max:
      return "\\/\\";
    case max_min:
      return "/\\/";
    case max:
      return "/\\";
    case max_left_min:
      return "/\\_";
    case right_min_max:
      return "_/\\";
  }
  AI_NEVER_REACHED
}

CubicEndShape get_end(CubicToNextSampleType type, bool left)
{
  using enum CubicToNextSampleType;
  switch (type)
  {
    case flat:                  // __
      return flat_low;
    case up:                    // /
      return uphill;
    case down:                  // \.
      return downhill;
    case right_stop:            // /^
      return left ? uphill : plus_inf;
    case left_stop:             // ^\.
      return left ? plus_inf : downhill;
    case right_min:             // _/
      return left ? flat_low : uphill;
    case left_min:              // \_
      return left ? downhill : flat_low;
    case right_max:             // ‾\.
      return left ? flat_high : downhill;
    case left_max:              // /‾
      return left ? uphill : flat_high;
    case right_max_left_min:    // ‾\_
      return left ? flat_high : flat_low;
    case right_min_left_max:    // _/‾
      return left ? flat_low : flat_high;
    case min:                   // \/
      return left ? downhill : uphill;
    case right_max_min:         // ‾\/
      return left ? flat_high : uphill;
    case min_left_max:          // \/‾
      return left ? downhill : flat_high;
    case min_max:               // \/\.
      return downhill;
    case max_min:               // /\/
      return uphill;
    case max:                   // /\.
      return left ? uphill : downhill;
    case max_left_min:          // /\_
      return left ? uphill : flat_low;
    case right_min_max:         // _/\.
      return left ? flat_low : downhill;
    default:
      ASSERT(false);
  }
  AI_NEVER_REACHED
}

} // namespace gradient_descent
