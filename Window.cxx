#include "sys.h"
#include "utils/AIAlert.h"
#include "Window.h"
#include "Layer.h"
#include "Plot.h"
#include "Draggable.h"
#include <X11/Xatom.h>
#include <mutex>
#include "debug.h"
#ifdef CWDEBUG
#include "debug_channel.h"
#include "debugcairo.h"
#endif

// Include this last as it defines a lot of macros that clash with C++.
#include <cairo/cairo-xlib.h>

#ifdef CWDEBUG
NAMESPACE_DEBUG_CHANNELS_START
channel_ct cairowindow("CAIROWINDOW");
NAMESPACE_DEBUG_CHANNELS_END
#endif

namespace cairowindow {

Window::Window(std::string title, int width, int height) : width_(width), height_(height), running_(false)
{
  DoutEntering(dc::cairowindow, "cairowindow::Window(\"" << title << "\", " << width << ", " << height << ") [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif

  display_ = XOpenDisplay(nullptr);
  if (!display_)
    THROW_ALERT("Cannot open DISPLAY.");

  int screen = DefaultScreen(display_);
  X11Window root = DefaultRootWindow(display_);

  // Create a window with title.
  x11window_ = XCreateSimpleWindow(display_, root, 0, 0, width, height, 1, BlackPixel(display_, screen), WhitePixel(display_, screen));
  XStoreName(display_, x11window_, title.c_str());

  // Register for WM_DELETE_WINDOW messages.
  wm_delete_window_ = XInternAtom(display_, "WM_DELETE_WINDOW", false);
  XSetWMProtocols(display_, x11window_, &wm_delete_window_, 1);

  custom_mouse_event_ = XInternAtom(display_, "CUSTOM_MOUSE_EVENT", false);

  // Create an X11 surface for the window.
  x11_surface_ = cairo_xlib_surface_create(display_, x11window_, DefaultVisual(display_, screen), width, height
    COMMA_CWDEBUG_ONLY("Window::x11_surface_:\"" + title + "\""));
  win_cr_ = cairo_create(x11_surface_
    COMMA_CWDEBUG_ONLY("Window::win_cr_:\"" + title + "\""));

  // Create an off-screen surface for double buffering.
  offscreen_surface_ = cairo_surface_create_similar(x11_surface_, CAIRO_CONTENT_COLOR, width, height
      COMMA_CWDEBUG_ONLY("Window::offscreen_surface_:\"" + title + "\""));
  offscreen_cr_ = cairo_create(offscreen_surface_
      COMMA_CWDEBUG_ONLY("Window::offscreen_cr_:\"" + title + "\""));
#ifdef CAIROWINDOW_DEBUGWINDOW
  debug_window_.start(offscreen_surface_, width, height, "offscreen_surface_");
#endif
}

Window::~Window()
{
  DoutEntering(dc::cairowindow, "Window::~Window() [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
#ifdef CAIROWINDOW_DEBUGWINDOW
  debug_window_.terminate();
#endif
  cairo_destroy(offscreen_cr_);
  cairo_surface_destroy(offscreen_surface_);
  cairo_destroy(win_cr_);
  cairo_surface_destroy(x11_surface_);
  if (mouse_grabbed_by_ != -1)
    release_mouse();
  XDestroyWindow(display_, x11window_);
  XCloseDisplay(display_);
}

void Window::add_printable(Printable* printable)
{
  DoutEntering(dc::notice, "Window::add_printable(" << (void*)printable << ") with geometry " << printable->geometry());
  printable_geometries_.emplace_back(printable->geometry(), printable);
}

// Called by XEventLoop thread.
void Window::grab_mouse(unsigned int button)
{
  DoutEntering(dc::notice, "Window::grab_mouse(" << button << ")");
  if ((mouse_button_mask_ & (1 << button)))
  {
    XGrabPointer(display_, x11window_, False,
        PointerMotionMask | ButtonReleaseMask, GrabModeAsync, GrabModeAsync, None, None, CurrentTime);
    mouse_grabbed_by_ = button;
  }
}

// Called by XEventLoop thread.
void Window::release_mouse()
{
  DoutEntering(dc::notice, "Window::release_mouse()");
  XUngrabPointer(display_, CurrentTime);
  mouse_grabbed_by_ = -1;
}

void Window::update(StrokeExtents const& stroke_extents)
{
  DoutEntering(dc::cairowindow, "Window::update(" << stroke_extents << ") [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif

  double const area_limit = stroke_extents.area() / 4;
  cairo_save(offscreen_cr_);
  stroke_extents.set_path(offscreen_cr_);
  cairo_clip(offscreen_cr_);
  {
    std::lock_guard<std::mutex> lock(offscreen_surface_mutex_);
    bool first_layer = true;
    cairo_set_operator(offscreen_cr_, CAIRO_OPERATOR_SOURCE);
    for (auto const& layer : layers_)
    {
      if (layer->area() < area_limit)
      {
        cairo_save(offscreen_cr_);
        layer->redraw(offscreen_cr_, stroke_extents);
        cairo_restore(offscreen_cr_);
      }
      else
      {
        cairo_set_source_surface(offscreen_cr_, layer->surface(), layer->offset_x(), layer->offset_y());
        cairo_paint(offscreen_cr_);
      }
      if (first_layer)
      {
        first_layer = false;
        cairo_set_operator(offscreen_cr_, CAIRO_OPERATOR_OVER);
      }
#ifdef CAIROWINDOW_DEBUGWINDOW
      debug_window_.update(stroke_extents);
#endif
    }
  }
  cairo_restore(offscreen_cr_);

  expose_events_.push_back(stroke_extents);

  if (send_expose_events_)
    set_send_expose_events(true);
}

void Window::set_send_expose_events(bool send_expose_events)
{
  send_expose_events_ = send_expose_events;

  if (send_expose_events_)
  {
    // Trigger an Expose event.
    XExposeEvent ev = {0};
    ev.type = Expose;
    ev.display = display_;
    ev.window = x11window_;

    for (int event = expose_events_.size() - 1; event >= 0; --event)
    {
      expose_events_[event].unpack(ev.x, ev.y, ev.width, ev.height);
      ev.count = event;
      XSendEvent(display_, x11window_, false, ExposureMask, (XEvent*)&ev);
    }
    expose_events_.clear();

    XFlush(display_);
  }
}

EventLoop Window::run()
{
  // Select input events.
  XSelectInput(display_, x11window_, ExposureMask | ButtonPressMask | ButtonReleaseMask | KeyPressMask | KeyReleaseMask);

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
  Debug(libcw_do.off());
  self->event_loop();
}

struct ExposeEventRect
{
  int x;
  int y;
  int width;
  int height;
};

void Window::event_loop()
{
  DoutEntering(dc::cairowindow, "Window::event_loop() [" << this << "]");
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  std::vector<ExposeEventRect> expose_events;
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
          if (!expose_events.empty())
          {
            for (ExposeEventRect const& rect : expose_events)
              cairo_rectangle(win_cr_, rect.x, rect.y, rect.width, rect.height);
            expose_events.clear();
          }
          cairo_rectangle(win_cr_, expose_event->x, expose_event->y, expose_event->width, expose_event->height);
          ASSERT(expose_event->x >= 0 && expose_event->y >= 0 && expose_event->x + expose_event->width <= 1200 && expose_event->y + expose_event->height <= 900);
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
      case KeyRelease:
      {
        XKeyEvent* key_event = (XKeyEvent*)&event;
        bool pressed = event.type == KeyPress;

        if (!pressed)
        {
          if (key_event->keycode == 107 || key_event->keycode == 65)
          {
            Message msg{InputEvent::key_release, key_event->x, key_event->y, key_event->keycode};
            push_message(msg);
          }
        }
        else
        {
          // Check if the PrtScn key was pressed.
          if (key_event->keycode == 107 || key_event->keycode == 65)
          {
            Message msg{InputEvent::key_press, key_event->x, key_event->y, key_event->keycode};
            push_message(msg);
            break;
          }
          running_ = false;           // Exit on any other key press.
          ++keypress_events;
        }
        break;
      }
      case ButtonPress:
      case ButtonRelease:
      {
        XButtonEvent* button_event = (XButtonEvent*)&event;
        bool pressed = event.type == ButtonPress;

        if (!pressed)
        {
          if (mouse_grabbed_by_ == event.xbutton.button)
            release_mouse();
          mouse_button_mask_ &= ~(1 << event.xbutton.button);
          Message msg{InputEvent::button_release, button_event->x, button_event->y, event.xbutton.button};
          push_message(msg);
        }
        else
        {
          mouse_button_mask_ |= (1 << event.xbutton.button);
          Message msg{InputEvent::button_press, button_event->x, button_event->y, event.xbutton.button};
          push_message(msg);
        }

        // Determine which button was pressed
        int button = button_event->button;
        Dout(dc::notice, "Button " << button << " was " << (pressed ? "pressed!" : "released!"));
        Dout(dc::notice, "x = " << button_event->x << "; y = " << button_event->y);
        Dout(dc::notice, "x_root = " << button_event->x_root << "; y_root = " << button_event->y_root);

        // Get the coordinates
        int x = button_event->x;
        int y = button_event->y;

        // Optionally, get root window coordinates
        int x_root = button_event->x_root;
        int y_root = button_event->y_root;

        break;
      }
      case MotionNotify:
      {
        if (mouse_grabbed_by_ != -1)
        {
          // Pass the mouse movement event to the main thread.
          Message msg{InputEvent::drag, event.xmotion.x, event.xmotion.y};
          push_message(msg);
        }
        break;
      }
      case ClientMessage:
      {
        if (static_cast<unsigned long>(event.xclient.data.l[0]) == wm_delete_window_)
          running_ = false;
        else if (static_cast<Atom>(event.xclient.message_type) == custom_mouse_event_)
        {
          Dout(dc::notice, "Received custom_mouse_event_");
          if (event.xclient.data.l[0] == custom_event_grab_mouse)
            grab_mouse(event.xclient.data.l[1]);
        }
        break;
      }
    }
  }
  Message msg{InputEvent::terminate_program};
  push_message(msg);
}

void Window::send_close_event()
{
  XEvent event = {};
  event.type = ClientMessage;
  event.xclient.window = x11window_;
  event.xclient.message_type = XInternAtom(display_, "WM_PROTOCOLS", true);
  event.xclient.format = 32;
  event.xclient.data.l[0] = wm_delete_window_;
  event.xclient.data.l[1] = CurrentTime;

  XSendEvent(display_, x11window_, false, NoEventMask, &event);
  XFlush(display_);
}

void Window::send_custom_event(uint32_t data, unsigned int button)
{
  XEvent event = {};
  event.xclient.type = ClientMessage;
  event.xclient.window = x11window_;
  event.xclient.message_type = custom_mouse_event_;
  event.xclient.format = 32;
  event.xclient.data.l[0] = data;
  event.xclient.data.l[1] = button;

  XSendEvent(display_, x11window_, false, NoEventMask, &event);
  XFlush(display_);
}

void Window::close()
{
  running_ = false;
  send_close_event();
}

void Window::register_draggable(plot::Plot& plot, plot::Draggable* draggable, std::function<Point (Point const&)> restriction)
{
  DoutEntering(dc::notice, "Window::register_draggable(@" << *draggable << ")");
  clickable_rectangles_.push_back(draggable->geometry());
  clickable_plots_.push_back(&plot);
  plot.register_draggable({}, draggable, std::move(restriction));
}

ClickableIndex Window::grab_draggable(double x, double y)
{
  DoutEntering(dc::notice, "Window::grab_draggable(" << x << ", " << y << ")");
  ClickableIndex found_index;
  double min_dist_squared = std::numeric_limits<double>::max();
  for (ClickableIndex index = clickable_rectangles_.ibegin(); index != clickable_rectangles_.iend(); ++index)
  {
    Rectangle const& geometry = clickable_rectangles_[index];
    // A Point uses ShapePosition::at_corner.
    double center_x = geometry.offset_x();
    double center_y = geometry.offset_y();
    double half_width = geometry.width();
    double half_height = geometry.height();
    if (center_x - half_width < x && x < center_x + half_width &&
        center_y - half_width < y && y < center_y + half_height)
    {
      double dist_squared = utils::square(center_x - x) + utils::square(center_y - y);
      if (dist_squared < min_dist_squared)
      {
        min_dist_squared = dist_squared;
        found_index = index;
      }
    }
  }
  return found_index;
}

bool Window::update_grabbed(ClickableIndex grabbed_point, double pixel_x, double pixel_y)
{
  plot::Plot* plot = clickable_plots_[grabbed_point];
  Rectangle new_rectangle = plot->update_grabbed({}, grabbed_point, pixel_x, pixel_y);
  if (new_rectangle.is_defined())
  {
    // Update the rectangle of a draggable Point, called after it was moved.
    clickable_rectangles_[grabbed_point] = new_rectangle;
    return true;
  }
  return false;
}

void Window::move_draggable(plot::Draggable* draggable, ClickableIndex clickable_index, cairowindow::Point new_position)
{
  plot::Plot* plot = clickable_plots_[clickable_index];
  draggable->set_position(new_position);        // Because we want to apply restrictions, if any, relative to the new position.
  plot->apply_restrictions({}, clickable_index, new_position);
  draggable->moved(plot, new_position);
  clickable_rectangles_[clickable_index] = draggable->geometry();
}

Printable* Window::find_printable(int mouse_x, int mouse_y)
{
  for (PrintableGeometries const& printable_geometries : printable_geometries_)
  {
    if (printable_geometries.geometry_.contains(mouse_x, mouse_y))
      return printable_geometries.printable_;
  }
  return nullptr;
}

bool Window::handle_input_events()
{
  bool block = true;
  while (have_message(block))
  {
    Message const* message = pop_message();

    Dout(dc::cairowindow, "Received message " << message->event << " (" << message->mouse_x << ", " << message->mouse_y << ")");
    switch (message->event)
    {
      case InputEvent::terminate_program:
      {
        return false;
      }
      case InputEvent::key_press:
      {
        Dout(dc::cairowindow, "key: " << message->detail.keycode);
        if (message->detail.keycode == 107)
        {
          Printable* printable = find_printable(message->mouse_x, message->mouse_y);
          if (printable)
            printable->set_need_print();
          return true;
        }
        else if (message->detail.keycode == 65)
        {
          return true;
        }
        break;
      }
      case InputEvent::key_release:
      {
        Dout(dc::cairowindow, "key: " << message->detail.keycode);
        break;
      }
      case InputEvent::button_press:
      {
        Dout(dc::cairowindow, "button: " << message->detail.button);
        auto index = grab_draggable(message->mouse_x, message->mouse_y);
        if (!index.undefined())
        {
          send_custom_event(custom_event_grab_mouse, message->detail.button);
          grab_index_ = index;
          grab_button_ = message->detail.button;
        }
        break;
      }
      case InputEvent::button_release:
      {
        Dout(dc::cairowindow, "button: " << message->detail.button);
        if (!grab_index_.undefined() && message->detail.button == grab_button_)
          grab_index_.set_to_undefined();
        break;
      }
      case InputEvent::drag:
        if (!grab_index_.undefined())
        {
          // Update the object with grab_index_ and return true if a redraw is necessary.
          if (update_grabbed(grab_index_, message->mouse_x, message->mouse_y))
            block = false;  // We have to redraw a part of the graph.
        }
        break;
    }
  }
  return true;
}

} // namespace cairowindow
