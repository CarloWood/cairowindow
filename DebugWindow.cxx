#include "sys.h"
#include "DebugWindow.h"
#include <chrono>
#include <iostream>
#include "debug.h"

namespace cairowindow {

cairo_t* DebugWindowVars::open_window(std::string const& title)
{
  display_ = XOpenDisplay(nullptr);
  if (!display_)
  {
    std::cerr << "Cannot open display.\n";
    return nullptr;
  }

  int screen = DefaultScreen(display_);
  Window root = DefaultRootWindow(display_);
  x11window_ = XCreateSimpleWindow(display_, root, 10, 10, width_, height_, 1, BlackPixel(display_, screen), WhitePixel(display_, screen));
  XStoreName(display_, x11window_, title.c_str());
  XMapWindow(display_, x11window_);
  XSelectInput(display_, x11window_, ExposureMask | KeyPressMask);

  x11_surface_ = cairo_xlib_surface_create(display_, x11window_, DefaultVisual(display_, screen), width_, height_);
  return cairo_create(x11_surface_);
}

void DebugWindowVars::close_window(cairo_t* cr)
{
  cairo_destroy(cr);
  cairo_surface_destroy(x11_surface_);
  XDestroyWindow(display_, x11window_);
  XCloseDisplay(display_);
}

void DebugWindow::main_loop(cairo_surface_t* shared_surface, std::string title)
{
  Debug(NAMESPACE_DEBUG::init_thread(title));

  cairo_t* cr = variables_type::wat{variables_}->open_window(title);
  until_window_opened_.open();

  running_ = true;
  while (running_)
  {
    for (;;)
    {
      XEvent event;
      {
        variables_type::wat vars_w(variables_);
        if (XPending(vars_w->display_) <= 0)
          break;
        XNextEvent(vars_w->display_, &event);
      }

      if (event.type == Expose)
      {
        std::lock_guard<std::mutex> lock(surface_mutex_);
        cairo_set_source_surface(cr, shared_surface, 0, 0);
        cairo_paint(cr);
      }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(200)); // Limit CPU usage.
  }

  variables_type::wat{variables_}->close_window(cr);
}

void DebugWindow::start(cairo_surface_t* shared_surface, double width, double height, std::string debug_name)
{
  {
    variables_type::wat variables_w(variables_);
    variables_w->width_ = width;
    variables_w->height_ = height;
  }

  debug_window_ = std::thread(thread_function, this, shared_surface, debug_name);
  until_window_opened_.wait();
}

void DebugWindow::terminate()
{
  running_ = false;
  debug_window_.join();
}

void DebugWindow::update(Rectangle const& rectangle)
{
  variables_type::wat variables_w(variables_);
  variables_w->trigger_expose_event();
}

void DebugWindowVars::trigger_expose_event()
{
  // Trigger an Expose event.
  XExposeEvent ev = {0};
  ev.type = Expose;
  ev.display = display_;
  ev.window = x11window_;
  ev.x = 0;
  ev.y = 0;
  ev.width = width_;
  ev.height = height_;
  ev.count = 0; // No more Expose events to follow.

  XSendEvent(display_, x11window_, false, ExposureMask, (XEvent*)&ev);
  XFlush(display_);
}

} // namespace cairowindow
