#include "sys.h"
#include "EventLoop.h"
#include "Window.h"
#include "debug.h"

namespace cairowindow {

EventLoop::EventLoop(Window* window, void (*thread_function)(Window*)) :
  cleanly_terminated_(false), window_(window), event_loop_thread_(thread_function, window)
{
}

EventLoop::~EventLoop()
{
  if (!cleanly_terminated_)
    window_->close();
  if (event_loop_thread_.joinable())
  {
    Dout(dc::notice, "Waiting for the window to be closed.");
    event_loop_thread_.join();
  }
}

} // namespace cairowindow
