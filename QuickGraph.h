#pragma once

#include "Range.h"
#include "Window.h"
#include "Plot.h"
#include "Line.h"
#include "LinePiece.h"
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
  Range y_range_;

  Window window_;
  boost::intrusive_ptr<Layer> background_layer_;
  boost::intrusive_ptr<Layer> second_layer_;
  std::thread event_loop_;
  plot::Plot plot_;
  std::vector<plot::BezierFitter> plot_bezier_fitter_;
  std::vector<plot::Line> plot_lines_;
  std::vector<plot::LinePiece> plot_line_pieces_;
  std::vector<plot::Point> plot_points_;

 public:
  QuickGraph(std::string const& title, std::string const& x_label, std::string const& y_label, Range x_range, Range y_range = {0.0, 0.0});

  QuickGraph(std::string const& title, std::string const& x_label, std::string const& y_label, Range x_range,
      std::function<double(double)> const& f, draw::LineStyle const& line_style) :
    QuickGraph(title, x_label, y_label, x_range, {0.0, 0.0})
  {
    add_function(f, line_style);
  }

  QuickGraph(std::string const& title, std::string const& x_label, std::string const& y_label, Range x_range, Range y_range,
      std::function<double(double)> const& f, draw::LineStyle const& line_style) :
    QuickGraph(title, x_label, y_label, x_range, y_range)
  {
    add_function(f, line_style);
  }

  QuickGraph(std::string const& title, std::string const& x_label, std::string const& y_label, Range x_range,
      std::function<double(double)> f, Color line_color = color::black) :
    QuickGraph(title, x_label, y_label, x_range, {0.0, 0.0})
  {
    add_function(f, line_color);
  }

  QuickGraph(std::string const& title, std::string const& x_label, std::string const& y_label, Range x_range, Range y_range,
      std::function<double(double)> f, Color line_color = color::black) :
    QuickGraph(title, x_label, y_label, x_range, y_range)
  {
    add_function(f, line_color);
  }

  ~QuickGraph();

  Range const& xrange() const { return x_range_; }
  Range const& yrange() const { return y_range_; }

  void add_function(std::function<double(double)> const& f, draw::LineStyle const& line_style);

  void add_function(std::function<double(double)> const& f, Color line_color = color::black)
  {
    draw::LineStyle line_style;
    add_function(f, line_style({.line_color = line_color, .line_width = 1.0}));
  }

  void add_line(Line const& L, draw::LineStyle const& line_style);

  void add_line(Line const& L, Color line_color = color::black)
  {
    draw::LineStyle line_style;
    add_line(L, line_style({.line_color = line_color, .line_width = 1.0}));
  }

  void add_line(LinePiece const& L, draw::LineStyle const& line_style);

  void add_line(LinePiece const& L, Color line_color = color::black)
  {
    draw::LineStyle line_style;
    add_line(L, line_style({.line_color = line_color, .line_width = 1.0}));
  }

  void add_point(Point P, draw::PointStyle const& point_style);

  void add_point(Point P)
  {
    draw::PointStyle point_style;
    add_point(P, point_style);
  }

  void clear()
  {
    plot_bezier_fitter_.clear();
    plot_lines_.clear();
    plot_points_.clear();
  }

  void wait_for_keypress();

 private:
  void initialize();
};

} // namespace cairowindow
