#pragma once

#include "AlgorithmEvent2.h"
#include "gradient_descent2/AlgorithmEventType.h"
#include "cairowindow/draw/Line.h"
#include "cairowindow/Window.h"
#include "cairowindow/Layer.h"
#include "cairowindow/Range.h"
#include "cairowindow/Plot.h"
#include "cairowindow/BezierFitter.h"
#include "events/Events.h"
#include <thread>
#include "debug.h"

namespace gradient_descent {
class Algorithm;
} // namespace gradient_descent

namespace enable_drawing {

class Function;

class EnableDrawing
{
 public:
  static constexpr cairowindow::draw::LineStyle curve_line_style{{.line_width = 1.0}};

 private:
  cairowindow::Window window;
  boost::intrusive_ptr<cairowindow::Layer> background_layer;
  boost::intrusive_ptr<cairowindow::Layer> second_layer_;
  std::thread event_loop;
  cairowindow::Range L_min_max;
  cairowindow::plot::Plot plot_;
  AlgorithmEvent algorithm_event;
  events::RequestHandle<gradient_descent::AlgorithmEventType> algorithm_event_handle;
  cairowindow::plot::BezierFitter plot_curve;

 public:
  EnableDrawing(gradient_descent::Algorithm* algorithm, Function const& L, double w_min, double w_max);
  ~EnableDrawing();

  cairowindow::plot::Plot& plot() { return plot_; }
  boost::intrusive_ptr<cairowindow::Layer> const& second_layer() const { return second_layer_; }

  static cairowindow::Range get_L_min_max(Function const& L, double w_min, double w_max);

  void wait()
  {
    // Block until a key is pressed.
    if (!window.handle_input_events())
      Dout(dc::notice, "Algorithm Terminated");
  }
};

} // namespace enable_drawing