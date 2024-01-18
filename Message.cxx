#include "sys.h"
#include "Message.h"
#include "utils/macros.h"
#include "debug.h"

namespace cairowindow {

#ifdef CWDEBUG
char const* to_string(MouseEvent mouse_event)
{
  switch (mouse_event)
  {
    AI_CASE_RETURN(MouseEvent::button_press);
    AI_CASE_RETURN(MouseEvent::button_release);
    AI_CASE_RETURN(MouseEvent::drag);
  }
  AI_NEVER_REACHED
}

std::ostream& operator<<(std::ostream& os, MouseEvent mouse_event)
{
  return os << to_string(mouse_event);
}

#endif

} // namespace cairowindow
