#pragma once

#include "Range.h"
#include "Window.h"
#include "Plot.h"
#include <string>
#include <functional>
#include <thread>

namespace cairowindow {

class QuickGraph
{
 private:
  bool empty_{true};
  std::string title_;
  std::string x_label_;
  std::string y_label_;
  Range x_range_;

  Window window_;
  boost::intrusive_ptr<Layer> background_layer_;
  boost::intrusive_ptr<Layer> second_layer_;
  std::thread event_loop_;
  plot::Plot plot_;
  std::vector<plot::BezierFitter> plot_bezier_fitter_;
  std::vector<plot::Point> plot_points_;

 public:
  QuickGraph(std::string const& title, std::string const& x_label, std::string const& y_label, Range x_range);

  QuickGraph(std::string const& title, std::string const& x_label, std::string const& y_label, Range x_range,
      std::function<double(double)> const& f, draw::LineStyle const& line_style) :
    QuickGraph(title, x_label, y_label, x_range)
  {
    add_function(f, line_style);
  }

  QuickGraph(std::string const& title, std::string const& x_label, std::string const& y_label, Range x_range,
      std::function<double(double)> f, Color line_color = color::black) :
    QuickGraph(title, x_label, y_label, x_range)
  {
    add_function(f, line_color);
  }

  ~QuickGraph();

  void add_function(std::function<double(double)> const& f, draw::LineStyle const& line_style);

  void add_function(std::function<double(double)> const& f, Color line_color = color::black)
  {
    draw::LineStyle line_style;
    add_function(f, line_style({.line_color = line_color, .line_width = 1.0}));
  }

  void add_point(Point P, draw::PointStyle const& point_style);

  void add_point(Point P)
  {
    draw::PointStyle point_style;
    add_point(P, point_style);
  }

  void wait_for_keypress();
};

} // namespace cairowindow
