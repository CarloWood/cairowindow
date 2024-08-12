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

} // namespace gradient_descent
