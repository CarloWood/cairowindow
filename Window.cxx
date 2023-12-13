#include "sys.h"
#include "utils/AIAlert.h"
#include "Window.h"
#include <X11/Xatom.h>
#include <mutex>
#include "debug.h"

namespace cairowindow {

Window::Window(std::string title, int width, int height) : width_(width), height_(height), running_(false)
{
  DoutEntering(dc::notice, "cairowindow::Window(\"" << title << "\", " << width << ", " << height << ")");

  display_ = XOpenDisplay(nullptr);
  if (!display_)
    THROW_ALERT("Cannot open DISPLAY.");

  int screen = DefaultScreen(display_);
  X11Window root = DefaultRootWindow(display_);

  // Create a window with title.
  x11window_ = XCreateSimpleWindow(display_, root, 0, 0, width, height, 1, BlackPixel(display_, screen), WhitePixel(display_, screen));
  XStoreName(display_, x11window_, title.c_str());

  // Register for WM_DELETE_WINDOW messages.
  wm_delete_window_ = XInternAtom(display_, "WM_DELETE_WINDOW", False);
  XSetWMProtocols(display_, x11window_, &wm_delete_window_, 1);

  // Create an X11 surface for the window.
  x11_surface_ = cairo_xlib_surface_create(display_, x11window_, DefaultVisual(display_, screen), width, height);
  win_cr_ = cairo_create(x11_surface_);

  // Create an off-screen surface for double buffering.
  offscreen_surface_ = cairo_surface_create_similar(x11_surface_, CAIRO_CONTENT_COLOR, width, height);
  offscreen_cr_ = cairo_create(offscreen_surface_);
}

Window::~Window()
{
  cairo_destroy(offscreen_cr_);
  cairo_surface_destroy(offscreen_surface_);
  cairo_destroy(win_cr_);
  cairo_surface_destroy(x11_surface_);
  XDestroyWindow(display_, x11window_);
  XCloseDisplay(display_);
}

Layer& Window::create_background_layer(Rectangle rectangle, Color background_color)
{
  if (!background_color.is_opaque())
    THROW_FALERT("The background layer can not have transparency.");
  layers_.emplace_back(x11_surface_, rectangle, CAIRO_CONTENT_COLOR, background_color, this);
  // Send a redraw event for the entire layer (because of the background color).
  redraw(rectangle);
  return layers_.back();
}

Layer& Window::create_layer(Rectangle rectangle)
{
  layers_.emplace_back(x11_surface_, rectangle, CAIRO_CONTENT_COLOR_ALPHA, Color{0, 0, 0, 0}, this);
  // Send a redraw event for the entire layer (because of the background color).
  redraw(rectangle);
  return layers_.back();
}

void Window::redraw(Rectangle const& rect)
{
  cairo_save(offscreen_cr_);
  cairo_rectangle(offscreen_cr_, rect.offset_x(), rect.offset_y(), rect.width(), rect.height());
  cairo_clip(offscreen_cr_);
  {
    std::lock_guard<std::mutex> lock(offscreen_surface_mutex_);
    for (Layer& layer : layers_)
    {
      cairo_set_source_surface(offscreen_cr_, layer.drawing_surface(), layer.offset_x(), layer.offset_y());
      cairo_paint(offscreen_cr_);
    }
  }
  cairo_restore(offscreen_cr_);

  // Trigger an Expose event.
  XExposeEvent ev = {0};
  ev.type = Expose;
  ev.display = display_;
  ev.window = x11window_;
  ev.x = rect.offset_x();
  ev.y = rect.offset_y();
  ev.width = rect.width();
  ev.height = rect.height();
  ev.count = 0; // No more Expose events to follow.

  XSendEvent(display_, x11window_, False, ExposureMask, (XEvent*)&ev);
  XFlush(display_);
}

EventLoop Window::run()
{
  // Select input events.
  XSelectInput(display_, x11window_, ExposureMask | ButtonPressMask | KeyPressMask);

  // Show the window.
  XMapWindow(display_, x11window_);

  // Start event loop thread.
  running_ = true;
  return {this, Window::event_loop_thread};
}

//static
void Window::event_loop_thread(Window* self)
{
  Debug(NAMESPACE_DEBUG::init_thread("XEventLoop"));
  self->event_loop();
}

struct ExposeEventRect
{
  int x;
  int y;
  int width;
  int height;
};

std::vector<ExposeEventRect> expose_events;

void Window::event_loop()
{
  int keypress_events = 0;
  // Event loop.
  XEvent event;
  while (running_)
  {
    // Block till the next X11 event.
    XNextEvent(display_, &event);
    switch (event.type)
    {
      case Expose:
      {
        XExposeEvent* expose_event = &event.xexpose;
        if (expose_event->count > 0)
          expose_events.emplace_back(expose_event->x, expose_event->y, expose_event->width, expose_event->height);
        else
        {
          cairo_save(win_cr_);
          if (expose_events.empty())
            cairo_rectangle(win_cr_, expose_event->x, expose_event->y, expose_event->width, expose_event->height);
          else
          {
            for (ExposeEventRect const& rect : expose_events)
              cairo_rectangle(win_cr_, rect.x, rect.y, rect.width, rect.height);
            expose_events.clear();
          }
          cairo_clip(win_cr_);
          {
            std::lock_guard<std::mutex> lock(offscreen_surface_mutex_);
            cairo_set_source_surface(win_cr_, offscreen_surface_, 0, 0);
            cairo_paint(win_cr_);
          }
          cairo_restore(win_cr_);
        }
        break;
      }
      case KeyPress:
      {
        running_ = false;       // Exit on any key press.
        ++keypress_events;
        break;
      }
      case ClientMessage:
      {
        if (static_cast<unsigned long>(event.xclient.data.l[0]) == wm_delete_window_)
          running_ = false;
        break;
      }
    }
  }
}

void Window::send_close_event()
{
  XEvent event = {};
  event.type = ClientMessage;
  event.xclient.window = x11window_;
  event.xclient.message_type = XInternAtom(display_, "WM_PROTOCOLS", True);
  event.xclient.format = 32;
  event.xclient.data.l[0] = wm_delete_window_;
  event.xclient.data.l[1] = CurrentTime;

  XSendEvent(display_, x11window_, False, NoEventMask, &event);
  XFlush(display_);
}

void Window::close()
{
  running_ = false;
  send_close_event();
}

} // namespace cairowindow
