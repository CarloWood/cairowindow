#include "sys.h"
#include "Message.h"
#include "utils/macros.h"
#include "debug.h"

namespace cairowindow {

#ifdef CWDEBUG
char const* to_string(InputEvent mouse_event)
{
  switch (mouse_event)
  {
    AI_CASE_RETURN(InputEvent::button_press);
    AI_CASE_RETURN(InputEvent::button_release);
    AI_CASE_RETURN(InputEvent::drag);
  }
  AI_NEVER_REACHED
}

std::ostream& operator<<(std::ostream& os, InputEvent mouse_event)
{
  return os << to_string(mouse_event);
}

#endif

} // namespace cairowindow
