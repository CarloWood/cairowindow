#pragma once

#include "Rectangle.h"
#include "StrokeExtents.h"
#include "EventLoop.h"
#include "LayerArgs.h"
#include "Message.h"
#include "Point.h"
#include "Draggable.h"
#include "plot/Point.h"
#include "utils/AIAlert.h"
#include "utils/threading/Semaphore.h"
#include "utils/threading/FIFOBuffer.h"
#include "utils/Vector.h"
#include <boost/intrusive_ptr.hpp>
#include <mutex>
#include <string>
#include <atomic>
#include <functional>
#include "debug.h"
#ifdef CAIROWINDOW_DEBUGWINDOW
#include "DebugWindow.h"
#endif
#ifdef CWDEBUG
#include "debug_channel.h"
#endif

// Get some X11 types that we need.
extern "C" {
#include <X11/Xdefs.h>                  // For Atom and XID
typedef XID X11Window;
typedef XID Drawable;
struct _XDisplay;                       // Forward declaration.
}

namespace cairowindow {
namespace plot {
class Plot;
class Draggable;
} // namespace plot

template<CS> class CoordinateSystem;

class Layer;
class Printable;

template<typename Type>
concept LayerType = std::is_base_of_v<Layer, Type>;

static const uint32_t custom_event_grab_mouse = 1;

class Window
{
 private:
  _XDisplay* display_;
  X11Window x11window_;
  int width_;
  int height_;

  cairo_surface_t* x11_surface_;
  cairo_t* win_cr_;
  cairo_surface_t* offscreen_surface_;
  cairo_t* offscreen_cr_;
  std::mutex offscreen_surface_mutex_;

  std::atomic_bool running_;
  std::atomic_bool destroyed_{false};
  Atom wm_delete_window_;
  Atom custom_mouse_event_;

  std::vector<boost::intrusive_ptr<Layer>> layers_;

  bool send_expose_events_{true};
  std::vector<StrokeExtents> expose_events_;

  // Mouse events (accessed by the XEventLoop thread).
  int mouse_grabbed_by_{-1};
  int mouse_button_mask_{0};
  utils::threading::Semaphore message_semaphore_{0};
  utils::threading::FIFOBuffer<1, Message> message_buffer_{64};

  // Dragging.
  utils::Vector<Geometry, ClickableIndex> clickable_rectangles_;
  utils::Vector<plot::Plot*, ClickableIndex> clickable_plots_;
  utils::Vector<std::function<Geometry(double, double)>, ClickableIndex> draggable_update_;     // Stores lambdas that convert the new position, (x,y) in pixels,
                                                                                                // to the new bounding box in pixels.
  ClickableIndex grab_index_;           // The index of the object (Point) that was grabbed.
  unsigned int grab_button_;            // The mouse button that did the grabbing (only valid when grab_index is not undefined).

  // Printing.
  struct PrintableGeometries
  {
    Geometry geometry_;
    Printable* printable_;
  };
  std::vector<PrintableGeometries> printable_geometries_;

#ifdef CAIROWINDOW_DEBUGWINDOW
  DebugWindow debug_window_;
#endif

 public:
  Window(std::string title, int width, int height);
  ~Window();

  void add_printable(Printable* printable);
  Printable* find_printable(int x, int y);

  EventLoop run();
  static void event_loop_thread(Window* self);
  void event_loop();

  void close();
  bool destroyed() const { return destroyed_; }

  Geometry geometry() const { return {0, 0, static_cast<double>(width_), static_cast<double>(height_)}; }

  template<LayerType LT, typename... ARGS>
  boost::intrusive_ptr<LT> create_layer(LayerArgs la, ARGS&&... args)
  {
    DoutEntering(dc::cairowindow, "Window::create_layer<" << libcwd::type_info_of<LT>().demangled_name() << ">(" << la <<
        join_more(", ", args...) << ") [" << this << "]");
    Geometry geometry = la.has_geometry() ? la.geometry() : this->geometry();
    boost::intrusive_ptr<LT> layer = new LT(x11_surface_, geometry, CAIRO_CONTENT_COLOR_ALPHA,
        Color{0, 0, 0, 0}, this, std::forward<ARGS>(args)...);
    layers_.push_back(layer);
    // Send an update request for the size of the entire layer (because of the background color).
    update(geometry);
    return layer;
  }

  template<LayerType LT>
  boost::intrusive_ptr<LT> create_layer()
  {
    return create_layer<LT>({});
  }

  template<LayerType LT, typename... ARGS>
  boost::intrusive_ptr<LT> create_background_layer(BackgroundLayerArgs la, ARGS&&... args)
  {
    DoutEntering(dc::cairowindow, "Window::create_background_layer<" << libcwd::type_info_of<LT>().demangled_name() << ">(" << la <<
        join_more(", ", args...) << ") [" << this << "]");
    if (!la.background_color().is_opaque())
      THROW_FALERT("The background layer can not have transparency.");
    Geometry geometry = la.has_geometry() ? la.geometry() : this->geometry();
    boost::intrusive_ptr<LT> layer = new LT(x11_surface_, geometry, CAIRO_CONTENT_COLOR,
        la.background_color(), this, std::forward<ARGS>(args)...);
    layer->add_area(geometry.area());
    layers_.push_back(layer);
    // Send an update request for the size of the entire layer (because of the background color).
    update(geometry);
    return layer;
  }

  // send_expose_events_ must be true for update to send expose events.
  void set_send_expose_events(bool send_expose_events);
  void update(StrokeExtents const& rectangle_list);

  bool push_message(Message const& message)
  {
    bool success = message_buffer_.push(&message);
    if (success)
      message_semaphore_.post();
    return !success;
  }

  bool have_message(bool block)
  {
    if (block)
    {
      message_semaphore_.wait();
      return true;
    }

    return message_semaphore_.try_wait();
  }

  Message const* pop_message()
  {
    return message_buffer_.pop();
  }

#ifdef CWDEBUG
  friend std::ostream& operator<<(std::ostream& os, Window const* window_ptr)
  {
    os << "Window*";
    return os;
  }
#endif

  void send_custom_event(uint32_t data, unsigned int button);

  // Allow dragging of point with the mouse.
  void register_draggable(plot::Plot& plot, plot::Draggable* draggable, std::function<Point (Point const&)> restriction = {});

  // Block until a point was dragged by the user, a key or button was pressed or released.
  // Returns true when a redraw is required, false when the program should be terminated.
  bool handle_input_events();

  // Called by Slider::set_value.
  bool update_grabbed(ClickableIndex grabbed_point, double pixel_x, double pixel_y);

  // Called by Point::move.
  void move_draggable(plot::Draggable* draggable, ClickableIndex clickable_index, Point new_position);

  template<CS cs>
  void register_draggable(CoordinateSystem<cs>& coordinate_system,
      plot::cs::Point<cs>* plot_point_cs,
      std::function<cs::Point<cs> (cs::Point<cs> const&)> restriction = {})
  {
    // Convert a new position (x,y), in pixels, to the relocated bounding box in pixels.
    auto update_grabbed_pixels = [&coordinate_system, plot_point_cs, restriction = std::move(restriction)](double pixel_x, double pixel_y) -> Geometry {
      return coordinate_system.update_grabbed(plot_point_cs, pixel_x, pixel_y, restriction);
    };
    register_draggable_impl(plot_point_cs, std::move(update_grabbed_pixels));
  }

 private:
  void send_close_event();
  void grab_mouse(unsigned int button);
  void release_mouse();

  ClickableIndex grab_draggable(double x, double y);

  void register_draggable_impl(plot::Draggable* draggable, std::function<Geometry(double, double)> update_pixels);
};

} // namespace cairowindow
