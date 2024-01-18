#pragma once

#include <iosfwd>

namespace cairowindow {

enum class MouseEvent
{
  button_press,
  button_release,
  drag
};

#ifdef CWDEBUG
char const* to_string(MouseEvent mouse_event);
std::ostream& operator<<(std::ostream& os, MouseEvent mouse_event);
#endif

struct Message
{
  MouseEvent event;
  int mouse_x;
  int mouse_y;
  unsigned int button;          // Only valid when event is button_press or button_release.
};

} // namespace cairowindow
