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
  cairo_t* cr;
  double width;
  double height;
  {
    variables_type::wat variables_w{variables_};
    cr = variables_w->open_window(title);
    width = variables_w->width_;
    height = variables_w->height_;
  }
  until_window_opened_.open();

  // Define the size of the tiles in the checkerboard pattern.
  int const tileSize = 20;
  double const lightGray = 0.8; // Light gray for one set of tiles.
  double const darkGray = 0.5;  // Dark gray for the other set of tiles.

  // Create a pattern surface to draw the checkerboard pattern.
  cairo_surface_t* pattern_surface = cairo_surface_create_similar(cairo_get_target(cr), CAIRO_CONTENT_COLOR, tileSize * 2, tileSize * 2);
  cairo_t* pattern_cr = cairo_create(pattern_surface);

  // Draw the checkerboard pattern.
  // Light gray squares.
  cairo_set_source_rgb(pattern_cr, lightGray, lightGray, lightGray);
  cairo_rectangle(pattern_cr, 0, 0, tileSize, tileSize);
  cairo_rectangle(pattern_cr, tileSize, tileSize, tileSize, tileSize);
  cairo_fill(pattern_cr);

  // Dark gray squares.
  cairo_set_source_rgb(pattern_cr, darkGray, darkGray, darkGray);
  cairo_rectangle(pattern_cr, tileSize, 0, tileSize, tileSize);
  cairo_rectangle(pattern_cr, 0, tileSize, tileSize, tileSize);
  cairo_fill(pattern_cr);

  // Set the pattern as the source for the main context.
  cairo_pattern_t* pattern = cairo_pattern_create_for_surface(pattern_surface);
  cairo_pattern_set_extend(pattern, CAIRO_EXTEND_REPEAT);

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
        cairo_set_source(cr, pattern);
        cairo_rectangle(cr, 0, 0, width, height);
        cairo_fill(cr);
        cairo_set_source_surface(cr, shared_surface, 0, 0);
        cairo_paint(cr);
      }
    }
    std::this_thread::sleep_for(std::chrono::milliseconds(200)); // Limit CPU usage.
  }

  // Clean up.
  cairo_pattern_destroy(pattern);
  cairo_destroy(pattern_cr);
  cairo_surface_destroy(pattern_surface);

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

void DebugWindow::update(StrokeExtents const& stroke_extents)
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
