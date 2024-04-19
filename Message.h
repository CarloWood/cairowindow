#pragma once

#include <iosfwd>

namespace cairowindow {

enum class InputEvent
{
  terminate_program,
  key_press,
  key_release,
  button_press,
  button_release,
  drag
};

#ifdef CWDEBUG
char const* to_string(InputEvent mouse_event);
std::ostream& operator<<(std::ostream& os, InputEvent mouse_event);
#endif

struct Message
{
  InputEvent event;
  int mouse_x;
  int mouse_y;
  unsigned int button;          // Only valid when event is button_press or button_release.
};

} // namespace cairowindow
