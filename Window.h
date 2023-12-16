#pragma once

#include "Rectangle.h"
#include "EventLoop.h"
#include "LayerArgs.h"
#include "utils/AIAlert.h"
#include <cairo/cairo-xlib.h>
#include <boost/intrusive_ptr.hpp>
#include <mutex>
#include <string>
#include <atomic>

#include <X11/Xlib.h>
#undef True
#undef False

namespace cairowindow {

using X11Window = ::Window;

class Layer;

template<typename Type>
concept LayerType = std::is_base_of_v<Layer, Type>;

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

  std::vector<boost::intrusive_ptr<Layer>> layers_;

 public:
  Window(std::string title, int width, int height);
  ~Window();

  EventLoop run();
  static void event_loop_thread(Window* self);
  void event_loop();

  void close();

  Rectangle get_rectangle() const { return {0, 0, static_cast<double>(width_), static_cast<double>(height_)}; }

  template<LayerType LT, typename... ARGS>
  boost::intrusive_ptr<LT> create_layer(LayerArgs la, ARGS&&... args)
  {
    Rectangle rectangle = la.has_rectangle() ? la.rectangle() : get_rectangle();
    boost::intrusive_ptr<LT> layer = new LT(x11_surface_, rectangle, CAIRO_CONTENT_COLOR_ALPHA,
        Color{0, 0, 0, 0}, this, std::forward<ARGS>(args)...);
    layers_.push_back(layer);
    // Send an update request for the size of the entire layer (because of the background color).
    update(rectangle);
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
    if (!la.background_color().is_opaque())
      THROW_FALERT("The background layer can not have transparency.");
    Rectangle rectangle = la.has_rectangle() ? la.rectangle() : get_rectangle();
    boost::intrusive_ptr<LT> layer = new LT(x11_surface_, rectangle, CAIRO_CONTENT_COLOR,
        la.background_color(), this, std::forward<ARGS>(args)...);
    layer->add_area(rectangle.area());
    layers_.push_back(layer);
    // Send an update request for the size of the entire layer (because of the background color).
    update(rectangle);
    return layer;
  }

  void update(Rectangle const& rectangle);

 private:
  void send_close_event();
};

} // namespace cairowindow
