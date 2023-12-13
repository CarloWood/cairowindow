#pragma once

#include "Layer.h"
#include "Rectangle.h"
#include "EventLoop.h"
#include <cairo/cairo-xlib.h>
#include <X11/Xlib.h>
#include <mutex>
#include <string>
#include <atomic>
#include <list>

namespace cairowindow {

using X11Window = ::Window;

class Layer;

class Window
{
 private:
  Display* display_;
  X11Window x11window_;
  int width_;
  int height_;

  cairo_surface_t* x11_surface_;
  cairo_t* win_cr_;
  cairo_surface_t* offscreen_surface_;
  cairo_t* offscreen_cr_;
  std::mutex offscreen_surface_mutex_;

  std::atomic_bool running_;
  Atom wm_delete_window_;

  std::list<Layer> layers_;

 public:
  Window(std::string title, int width, int height);
  ~Window();

  EventLoop run();
  static void event_loop_thread(Window* self);
  void event_loop();

  void close();

  Rectangle get_rect() const { return {0, 0, static_cast<double>(width_), static_cast<double>(height_)}; }
  Layer& create_background_layer(Rectangle rectangle, Color background_color);
  Layer& create_background_layer(Color background_color) { return create_background_layer(get_rect(), background_color); }
  Layer& create_layer(Rectangle rectangle);
  Layer& create_layer() { return create_layer(get_rect()); }

  void redraw(Rectangle const& rect);

 private:
  void send_close_event();
};

} // namespace cairowindow
