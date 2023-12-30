#pragma once

#ifndef CAIROWINDOW_DEBUGWINDOW
#error "Do not include DebugWindow.h without setting cmake option -DEnableDebugWindow:BOOL=ON"
#endif

#include "StrokeExtents.h"
#include "threadsafe/threadsafe.h"
#include "utils/threading/Gate.h"
#include <cairo/cairo.h>
#include <mutex>
#include <string>
#include <atomic>
#include <thread>

#include <cairo/cairo-xlib.h>
#undef True
#undef False
#undef Status

namespace cairowindow {

class Rectangle;

struct DebugWindowVars
{
  using Window = ::Window;

  Display* display_;
  Window x11window_;
  cairo_surface_t* x11_surface_;
  double width_;
  double height_;

  [[nodiscard]] cairo_t* open_window(std::string const& title);
  void close_window(cairo_t* cr);

  void trigger_expose_event();
};

struct DebugWindow
{
  std::mutex surface_mutex_;    // Not used.
  std::thread debug_window_;
  std::atomic_bool running_;
  utils::threading::Gate until_window_opened_;

  using variables_type = threadsafe::Unlocked<DebugWindowVars, threadsafe::policy::Primitive<std::mutex>>;
  variables_type variables_;

  void main_loop(cairo_surface_t* shared_surface, std::string title);

  static void thread_function(DebugWindow* self, cairo_surface_t* shared_surface, std::string title)
  {
    Debug(NAMESPACE_DEBUG::init_thread(title));
    self->main_loop(shared_surface, title);
  }

  void start(cairo_surface_t* shared_surface, double width, double height, std::string debug_name);
  void terminate();

  void update(StrokeExtents const& stroke_extents);
};

} // namespace cairowindow
