#pragma once

#include "Rectangle.h"
#include "StrokeExtents.h"
#include "EventLoop.h"
#include "LayerArgs.h"
#include "Message.h"
#include "Point.h"
#include "Draggable.h"
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
using ClickableIndex = utils::VectorIndex<Rectangle>;
namespace plot {
class Plot;
namespace cs {
template<CS> class Point;
} // namespace cs
} // namespace plot

template<CS cs> class CoordinateMapper;

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
  utils::Vector<std::function<Geometry(math::cs::Point<csid::pixels>)>, ClickableIndex> draggable_update_;      // Stores lambdas that convert the new position, (x,y) in pixels,
                                                                                                                // to the new bounding box in pixels.
  utils::Vector<std::function<void(ClickableIndex, double, double)>, ClickableIndex> draggable_update_cs_;      // Stores lambdas that convert the new position, (x,y) in cs, to
                                                                                                                // pixels and calls update_grabbed with that.
  utils::Vector<Geometry, ClickableIndex> clickable_geometries_;

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

  // Block until a point was dragged by the user, a key or button was pressed or released.
  // Returns true when a redraw is required, false when the program should be terminated.
  bool handle_input_events();

  // Called by Slider::set_value.
  void update_grabbed(ClickableIndex clickable_index, math::cs::Point<csid::pixels> new_position_pixels);

  // Called by Point::move_to.
  template<CS cs>
  void move_draggable(plot::cs::Point<cs> const* draggable, ClickableIndex clickable_index, math::cs::Point<cs> new_position_cs)
  {
    // Fake dragging by calling update_grabbed.
    draggable_update_cs_[clickable_index](clickable_index, new_position_cs.x(), new_position_cs.y());
  }

  // Update the clickable geometry for a draggable after it was moved programmatically.
  // This is needed when the draggable did not move through Window::update_grabbed (dragging) or Window::move_draggable (plot::Point::move_to).
  void update_draggable_geometry(plot::Draggable const& draggable);

  template<CS cs>
  using RestrictionFunction = std::function<math::cs::Point<cs> (math::cs::Point<cs> const&)>;

  // Allow dragging of a point with the mouse.
  // Note: neither coordinate_mapper nor plot_point_cs, passed from the caller, maybe be destructed or moved in memory for as long as this draggable is active.
  template<CS cs, typename L = std::nullptr_t>
  void register_draggable(CoordinateMapper<cs>& coordinate_mapper,
      plot::cs::Point<cs>* plot_point_cs,
      L const& restriction_cs = nullptr)
  requires (
    std::is_null_pointer_v<std::remove_cvref_t<L>> ||
    requires { { std::declval<L const&>() } -> std::convertible_to<RestrictionFunction<cs>>; }
  )
  {
    RestrictionFunction<cs> restriction;
    if constexpr (!std::is_same_v<L, std::nullptr_t>)
      restriction = static_cast<RestrictionFunction<cs>>(restriction_cs);
    // Convert a new position (x,y) in pixels, to the relocated bounding box in pixels, by calling update_grabbed_cs on the coordinate_mapper.
    // CoordinateMapper<cs>::update_grabbed_cs converts the position to cs before passing it to restriction_cs and storing the new position in plot_point_cs.
    auto update_grabbed_pixels = [&coordinate_mapper, plot_point_cs, restriction = std::move(restriction)](math::cs::Point<csid::pixels> position_pixels) -> Geometry {
      return coordinate_mapper.update_grabbed_cs(plot_point_cs, position_pixels, restriction);
    };
    // Convert a moved to position (x,y) in cs to pixels and ...
    auto update_cs = [this, &coordinate_mapper](ClickableIndex index, double position_x_cs, double position_y_cs) {
      math::cs::Point<cs> const position_cs{position_x_cs, position_y_cs};
      math::cs::Point<csid::pixels> new_position_pixels = position_cs * coordinate_mapper.cs_transform_pixels();
      update_grabbed(index, new_position_pixels);
    };
    // Store the lambdas and clickable geometry in a Vector.
    ClickableIndex index = register_draggable_impl(std::move(update_grabbed_pixels), std::move(update_cs), plot_point_cs->geometry());
    plot_point_cs->set_index(index);
  }

  // Allow dragging of a draggable object that already lives in pixel coordinate space (for example draw::Slider).
  template<typename L = std::nullptr_t>
  void register_draggable(plot::Draggable* draggable,
      L const& restriction_pixels = nullptr)
  requires (
    std::is_null_pointer_v<std::remove_cvref_t<L>> ||
    requires { { std::declval<L const&>() } -> std::convertible_to<RestrictionFunction<csid::pixels>>; }
  )
  {
    RestrictionFunction<csid::pixels> restriction;
    if constexpr (!std::is_same_v<L, std::nullptr_t>)
      restriction = static_cast<RestrictionFunction<csid::pixels>>(restriction_pixels);

    // Convert a new position (x,y) in pixels to the relocated bounding box in pixels by calling Draggable::moved.
    auto update_grabbed_pixels = [draggable, restriction = std::move(restriction)](math::cs::Point<csid::pixels> position_pixels) -> Geometry {
      if (restriction)
        position_pixels = restriction(position_pixels);
      draggable->moved(position_pixels);
      return draggable->geometry();
    };
    // Convert a moved to position (x,y) in cs to pixels... but since cs == csid::pixels, that just comes down to calling update_grabbed.
    auto update_cs = [this](ClickableIndex index, double position_x_pixels, double position_y_pixels) {
      math::cs::Point<csid::pixels> new_position_pixels{position_x_pixels, position_y_pixels};
      update_grabbed(index, new_position_pixels);
    };
    // Store the lambdas and clickable geometry in a Vector.
    ClickableIndex index = register_draggable_impl(std::move(update_grabbed_pixels), std::move(update_cs), draggable->geometry());
    draggable->set_index(index);
  }

 private:
  void send_close_event();
  void grab_mouse(unsigned int button);
  void release_mouse();

  ClickableIndex grab_draggable(math::cs::Point<csid::pixels> const& mouse_position);
  ClickableIndex register_draggable_impl(
      std::function<Geometry (math::cs::Point<csid::pixels>)> update_grabbed_pixels,
      std::function<void (ClickableIndex index, double x_cs, double y_cs)> update_cs,
      Geometry const& current_geometry);
};

} // namespace cairowindow
