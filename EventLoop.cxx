#include "sys.h"
#include "EventLoop.h"
#include "Window.h"
#include "debug.h"

namespace cairowindow {

EventLoop::EventLoop(Window* window, void (*thread_function)(Window*)) :
  cleanly_terminated_(false), window_(window), event_loop_thread_(thread_function, window)
{
  DoutEntering(dc::cairowindow, "EventLoop::EventLoop(" << window << ", <event_loop_thread_ lambda>)");
}

EventLoop::~EventLoop()
{
  DoutEntering(dc::cairowindow, "EventLoop::~EventLoop()");
  if (!cleanly_terminated_)
    window_->close();
  if (event_loop_thread_.joinable())
  {
#ifdef CWDEBUG
    auto id = event_loop_thread_.get_id();
#endif
    Dout(dc::cairowindow, "Joining with thread " << std::hex << id << ": waiting for the window to be closed.");
    event_loop_thread_.join();
    Dout(dc::cairowindow, "Successfully joined with event_loop_thread_ (" << std::hex << id << ")!");
  }
}

} // namespace cairowindow
