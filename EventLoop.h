#pragma once

#include <thread>

namespace cairowindow {

class Window;

// Returned by Window::run().
// Call set_cleanly_terminated() before leaving its scope (other than by exception).
class EventLoop
{
 private:
  bool cleanly_terminated_;
  Window* window_;
  std::thread event_loop_thread_;

 public:
  EventLoop(Window* window, void (*thread_function)(Window*));
  EventLoop(EventLoop&& event_loop) = default;
  ~EventLoop();

  void set_cleanly_terminated() { cleanly_terminated_ = true; }
  bool cleanly_terminated() const { return cleanly_terminated_; }
};

} // namespace cairowindow
