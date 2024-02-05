#pragma once

#include "draw/PlotArea.h"
#include "draw/Text.h"
#include "draw/Point.h"
#include "draw/Circle.h"
#include "draw/Curve.h"
#include "draw/Connector.h"
#include "draw/Rectangle.h"
#include "draw/Arc.h"
#include "draw/BezierCurve.h"
#include "draw/Slider.h"
#include "Range.h"
#include "Point.h"
#include "Circle.h"
#include "Text.h"
#include "LinePiece.h"
#include "Line.h"
#include "Curve.h"
#include "Connector.h"
#include "Rectangle.h"
#include "Arc.h"
#include "BezierCurve.h"
#include "Draggable.h"
#include "utils/Vector.h"
#include "utils/Badge.h"
#include <boost/intrusive_ptr.hpp>
#include <string>
#include <vector>
#include <functional>
#include <algorithm>
#include <type_traits>
#include <tuple>
#include <utility>
#include "debug.h"
#ifdef CWDEBUG
#include "debug_channel.h"
#endif

namespace cairowindow {
class Window;

namespace plot {
using utils::has_print_on::operator<<;

class Plot;

//--------------------------------------------------------------------------
// Point

class Point : public cairowindow::Point, public Draggable
{
 public:
  Point() = default;
  Point(cairowindow::Point const& point) : cairowindow::Point(point) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Point> draw_object_;

  // Implementation of Draggable.
  cairowindow::Rectangle const& geometry() const override { return draw_object_->geometry(); }
  void moved(Plot* plot, cairowindow::Point const& new_position) override;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const override { cairowindow::Point::print_on(os); }
#endif
};

//--------------------------------------------------------------------------
// LinePiece

class LinePiece : public cairowindow::LinePiece
{
 public:
  using cairowindow::LinePiece::LinePiece;
  LinePiece(cairowindow::LinePiece const& line_piece) : cairowindow::LinePiece(line_piece) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

//--------------------------------------------------------------------------
// Connector

class Connector : public cairowindow::Connector
{
 public:
  using cairowindow::Connector::Connector;
  Connector(cairowindow::Connector const& connector) : cairowindow::Connector(connector) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Connector> draw_object_;
};

//--------------------------------------------------------------------------
// Line

class Line : public cairowindow::Line
{
 public:
  using cairowindow::Line::Line;
  Line(cairowindow::Line const& line) : cairowindow::Line(line) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

//--------------------------------------------------------------------------
// Rectangle

class Rectangle : public cairowindow::Rectangle
{
 public:
  using cairowindow::Rectangle::Rectangle;
  Rectangle(cairowindow::Rectangle const& rectangle) : cairowindow::Rectangle(rectangle) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Rectangle> draw_object_;
};

//--------------------------------------------------------------------------
// Circle

class Circle : public cairowindow::Circle
{
 public:
  using cairowindow::Circle::Circle;
  Circle(cairowindow::Circle const& circle) : cairowindow::Circle(circle) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Circle> draw_object_;
};

//--------------------------------------------------------------------------
// Arc

class Arc : public cairowindow::Arc
{
 public:
  using cairowindow::Arc::Arc;
  Arc(cairowindow::Arc const& arc) : cairowindow::Arc(arc) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Arc> draw_object_;
};

//--------------------------------------------------------------------------
// BezierCurve

class BezierCurve : public cairowindow::BezierCurve
{
 public:
  using cairowindow::BezierCurve::BezierCurve;
  BezierCurve(cairowindow::BezierCurve const& bezier_curve) : cairowindow::BezierCurve(bezier_curve) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::BezierCurve> draw_object_;
};

//--------------------------------------------------------------------------
// Text

class Text : public cairowindow::Text
{
 public:
  using cairowindow::Text::Text;
  Text(cairowindow::Text const& text) : cairowindow::Text(text) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Text> draw_object_;
};

//--------------------------------------------------------------------------
// Slider

class Slider
{
 private:
  cairowindow::Rectangle geometry_;
  double min_value_;      // The value corresponding to position 0.
  double max_value_;      // The value corresponding to position 1.

 public:
  Slider(cairowindow::Rectangle geometry, double min_value, double max_value) :
    geometry_(geometry), min_value_(min_value), max_value_(max_value) { }

  cairowindow::Rectangle const& geometry() const { return geometry_; }
  double value() const { return min_value_ + draw_object_->rel_value() * (max_value_ - min_value_); }
  double min_value() const { return min_value_; }
  double max_value() const { return max_value_; }

  void set_value(double value);

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Slider> draw_object_;
};

//--------------------------------------------------------------------------

class Curve : public cairowindow::Curve
{
 public:
  using cairowindow::Curve::Curve;
  Curve(std::vector<cairowindow::Point> const& points, std::shared_ptr<draw::Curve> const& draw_object) :
    cairowindow::Curve(points), draw_object_(draw_object) { }
  Curve(std::vector<cairowindow::Point>&& points, std::shared_ptr<draw::Curve> const& draw_object) :
    cairowindow::Curve(std::move(points)), draw_object_(draw_object) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Curve> draw_object_;
};

enum class LineExtend
{
  none = 0,
  from = 1,
  to = 2,
  both = from|to
};

namespace {

template<typename Tuple, std::size_t... Indices>
auto tuple_tail_impl(Tuple&& tuple, std::index_sequence<Indices...>)
{
  return std::make_tuple(std::get<Indices + 1>(std::forward<Tuple>(tuple))...);
}

template<typename... Args>
auto tuple_tail(std::tuple<Args...>&& tuple)
{
  return tuple_tail_impl(std::forward<std::tuple<Args...>>(tuple), std::make_index_sequence<sizeof...(Args) - 1>{});
}

} // namespace

class Plot
{
 private:
  static constexpr int number_of_axes = draw::PlotArea::number_of_axes;

  draw::PlotArea plot_area_;
  std::shared_ptr<draw::Text> title_;
  std::shared_ptr<draw::Text> xlabel_;
  std::shared_ptr<draw::Text> ylabel_;
  std::array<Range, number_of_axes> range_;
  std::array<int, number_of_axes> range_ticks_{{10, 10}};
  std::array<std::vector<std::shared_ptr<draw::Text>>, number_of_axes> labels_;

  utils::Vector<Draggable*, ClickableIndex> draggables_;
  utils::Vector<std::function<cairowindow::Point (cairowindow::Point const&)>, ClickableIndex> draggable_restrictions_;

  struct TitleStyleDefaults : draw::DefaultTextStyleDefaults
  {
    static constexpr draw::TextPosition position = draw::centered;
    static constexpr double font_size = 24.0;
    static constexpr double offset = 0.0;
  };

  struct LabelStyleDefaults : draw::DefaultTextStyleDefaults
  {
    static constexpr draw::TextPosition position = draw::centered_below;
    static constexpr double font_size = 18.0;
  };

  struct XLabelStyleDefaults : LabelStyleDefaults
  {
    static constexpr double offset = 10.0;
  };

  struct YLabelStyleDefaults : LabelStyleDefaults
  {
    static constexpr draw::TextPosition position = draw::centered_above;
    static constexpr double rotation = -0.5 * M_PI;
    static constexpr double offset = XLabelStyleDefaults::offset;
  };

 public:
  using TitleStyle = draw::TextStyle<TitleStyleDefaults>;
  using XLabelStyle = draw::TextStyle<XLabelStyleDefaults>;
  using YLabelStyle = draw::TextStyle<YLabelStyleDefaults>;
  using LabelStyle = draw::TextStyle<LabelStyleDefaults>;

 public:
  Plot(cairowindow::Rectangle const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, TitleStyle title_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width), plot_area_style),
    title_(std::make_shared<draw::Text>(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_style.offset, title_style)) { }

  Plot(cairowindow::Rectangle const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, TitleStyle title_style,
      std::string xlabel, XLabelStyle xlabel_style, std::string ylabel, YLabelStyle ylabel_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width), plot_area_style),
    title_(std::make_shared<draw::Text>(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_style.offset, title_style)),
    xlabel_(std::make_shared<draw::Text>(xlabel, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() + plot_area_.geometry().height() + XLabelStyleDefaults::offset, xlabel_style)),
    ylabel_(std::make_shared<draw::Text>(ylabel, plot_area_.geometry().offset_x() - YLabelStyleDefaults::offset,
        plot_area_.geometry().offset_y() + 0.5 * plot_area_.geometry().height(), ylabel_style)) { }

  Plot(cairowindow::Rectangle const& geometry, draw::PlotAreaStyle plot_area_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width), plot_area_style) { }

  void set_range(int axis, Range range)
  {
    DoutEntering(dc::cairowindow, "Plot::set_range(" << axis << ", " << range << ") [" << this << "]");
    range_[axis] = range;
    range_ticks_[axis] = draw::PlotArea::calculate_range_ticks(range_[axis]);
    Dout(dc::cairowindow, "range_[" << axis << "] = " << range_[axis] << "; range_ticks_[" << axis << "] = " << range_ticks_[axis]);
  }

  void set_xrange(Range x_range) { set_range(x_axis, x_range); }
  void set_yrange(Range y_range) { set_range(y_axis, y_range); }

  cairowindow::Point clamp_to_plot_area(cairowindow::Point const& point) const
  {
    return {std::clamp(point.x(), range_[x_axis].min(), range_[x_axis].max()),
            std::clamp(point.y(), range_[y_axis].min(), range_[y_axis].max())};
  }

  double convert_x(double x) const;
  double convert_y(double y) const;
  Pixel convert_to_pixel(cairowindow::Point const& point) const;

  //--------------------------------------------------------------------------
  // Point

  // Add and draw plot_point on layer using point_style.
  void add_point(boost::intrusive_ptr<Layer> const& layer,
      draw::PointStyle const& point_style,
      Point const& plot_point);

  // Create and draw a point on layer at x,y using point_style.
  [[nodiscard]] Point create_point(boost::intrusive_ptr<Layer> const& layer,
      draw::PointStyle const& point_style,
      cairowindow::Point const& point)
  {
    Point plot_point(point);
    add_point(layer, point_style, plot_point);
    return plot_point;
  }

  //--------------------------------------------------------------------------
  // LinePiece

  // Add and draw plot_line_piece on layer using line_style and line_extend.
  void add_line(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style, LineExtend line_extend,
      LinePiece const& plot_line_piece);

  // Create and draw a line piece between points from and to using line_style and line_extend.
  template<typename... Args>
  [[nodiscard]] LinePiece create_line(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style, LineExtend line_extend,
      Args&&... args)
    requires requires(Args&&... args) { LinePiece{std::forward<Args>(args)...}; }
  {
    LinePiece plot_line_piece(std::forward<Args>(args)...);
    add_line(layer, line_style, line_extend, plot_line_piece);
    return plot_line_piece;
  }

  // Same, but without a line_extend.
  template<typename... Args>
  [[nodiscard]] LinePiece create_line(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      Args&&... args)
    requires requires(Args&&... args) { LinePiece{std::forward<Args>(args)...}; }
  {
    LinePiece plot_line_piece(std::forward<Args>(args)...);
    add_line(layer, line_style, LineExtend::none, plot_line_piece);
    return plot_line_piece;
  }

  //--------------------------------------------------------------------------
  // Connector

  // Add and draw plot_connector on layer, using line_style, fill_color, arrow_head_shape_from and arrow_head_shape_to.
  void add_connector(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style, Color fill_color,
      Connector::ArrowHeadShape arrow_head_shape_from, Connector::ArrowHeadShape arrow_head_shape_to,
      Connector const& plot_connector);

 private:
  template<typename... Args>
  [[nodiscard]] Connector create_connector_helper(
      Color& fill_color,
      Connector::ArrowHeadShape& arrow_head_shape_from,
      Connector::ArrowHeadShape& arrow_head_shape_to,
      boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style, std::tuple<Args...>&& args)
  {
    using type_of_first_arg = std::tuple_element_t<0, std::tuple<Args...>>;
    if constexpr (std::is_same_v<type_of_first_arg, Color> ||
                  std::is_same_v<type_of_first_arg, typename Connector::ArrowHeadShape>)
    {
      if constexpr (std::is_same_v<type_of_first_arg, Color>)
        fill_color = std::get<0>(args);
      else
      {
        if constexpr (std::is_same_v<std::tuple_element_t<1, std::tuple<Args...>>, typename Connector::ArrowHeadShape>)
          arrow_head_shape_from = std::get<0>(args);
        else
          arrow_head_shape_to = std::get<0>(args);
      }
      return create_connector_helper(fill_color, arrow_head_shape_from, arrow_head_shape_to, layer, line_style, tuple_tail(std::move(args)));
    }
    else
    {
      Connector plot_connector = std::apply([&](auto&&... unpacked_args) -> Connector {
        return {std::forward<decltype(unpacked_args)>(unpacked_args)..., arrow_head_shape_from, arrow_head_shape_to};
      }, std::move(args));

      add_connector(layer, line_style, fill_color, arrow_head_shape_from, arrow_head_shape_to, plot_connector);
      return plot_connector;
    }
  }

 public:
  // Create and draw a connector on layer, using args... and line_style.
  // Args can optionally start with a fill_color and/or zero, one or two ArrowHeadShape arguments.
  template<typename... Args>
  [[nodiscard]] Connector create_connector(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style, Args&&... args)
  {
    Color fill_color = color::white;
    Connector::ArrowHeadShape arrow_head_shape_from = Connector::no_arrow;
    Connector::ArrowHeadShape arrow_head_shape_to = Connector::open_arrow;
    return create_connector_helper(fill_color, arrow_head_shape_from, arrow_head_shape_to, layer,
        line_style, std::make_tuple(std::forward<Args>(args)...));
  }

  //--------------------------------------------------------------------------
  // Line

  // Add and draw plot_line using line_style.
  void add_line(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      Line const& plot_line);

  // Create and draw a line through point in direction using line_style.
  template<typename... Args>
  [[nodiscard]] Line create_line(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      Args&&... args)
  {
    Line plot_line(std::forward<Args>(args)...);
    add_line(layer, line_style, plot_line);
    return plot_line;
  }

  //--------------------------------------------------------------------------
  // Rectangle

  // Add and draw plot_rectangle on layer, using rectangle_style.
  void add_rectangle(boost::intrusive_ptr<Layer> const& layer, draw::RectangleStyle const& rectangle_style, Rectangle const& plot_rectangle);

  // Create and draw a rectangle on layer, using args... and rectangle_style.
  template<typename... Args>
  [[nodiscard]] Rectangle create_rectangle(boost::intrusive_ptr<Layer> const& layer,
      draw::RectangleStyle const& rectangle_style, Args&&... args)
  {
    Rectangle plot_rectangle(args...);
    add_rectangle(layer, rectangle_style, plot_rectangle);
    return plot_rectangle;
  }

  //--------------------------------------------------------------------------
  // Circle

  // Add and draw plot_circle on layer with center and radius using circle_style.
  void add_circle(boost::intrusive_ptr<Layer> const& layer,
      draw::CircleStyle const& circle_style,
      Circle const& plot_circle);

  // Create and draw a circle on layer with center and radius using circle_style.
  template<typename... Args>
  [[nodiscard]] Circle create_circle(boost::intrusive_ptr<Layer> const& layer,
      draw::CircleStyle const& circle_style,
      Args&&... args)
  {
    Circle plot_circle(std::forward<Args>(args)...);
    add_circle(layer, circle_style, plot_circle);
    return plot_circle;
  }

  // Same as above but use line_style (no fill_color).
  template<typename... Args>
  [[nodiscard]] Circle create_circle(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      Args&&... args)
  {
    return create_circle(layer, draw::CircleStyle({.line_color = line_style.line_color(), .line_width = line_style.line_width()}),
        std::forward<Args>(args)...);
  }

  //--------------------------------------------------------------------------
  // Arc

  // Add and draw plot_arc on layer using arc_style.
  void add_arc(boost::intrusive_ptr<Layer> const& layer,
      draw::ArcStyle const& arc_style,
      Arc const& plot_arc);

  // Create and draw an arc on layer width center, radius and start- and end_angle, using arc_style.
  template<typename... Args>
  [[nodiscard]] Arc create_arc(boost::intrusive_ptr<Layer> const& layer,
      draw::ArcStyle const& arc_style,
      Args&&... args)
  {
    Arc plot_arc(std::forward<Args>(args)...);
    add_arc(layer, arc_style, plot_arc);
    return plot_arc;
  }

  //--------------------------------------------------------------------------
  // BezierCurve

  // Add and draw plot_bezier_curve on layer using bezier_style.
  void add_bezier_curve(boost::intrusive_ptr<Layer> const& layer,
      draw::BezierCurveStyle const& bezier_style,
      BezierCurve const& plot_bezier_curve);

  // Create and draw a Bezier curve on layer using bezier_style.
  template<typename... Args>
  [[nodiscard]] BezierCurve create_bezier_curve(boost::intrusive_ptr<Layer> const& layer,
      draw::BezierCurveStyle const& bezier_style,
      Args&&... args)
  {
    BezierCurve plot_bezier_curve(std::forward<Args>(args)...);
    add_bezier_curve(layer, bezier_style, plot_bezier_curve);
    return plot_bezier_curve;
  }

  //--------------------------------------------------------------------------
  // Text

  // Add and draw plot_text on layer using text_style.
  void add_text(boost::intrusive_ptr<Layer> const& layer,
      draw::TextStyle<> const& text_style,
      Text const& plot_text);

  // Create and draw text on layer at position using text_style.
  [[nodiscard]] Text create_text(boost::intrusive_ptr<Layer> const& layer,
      draw::TextStyle<> const& text_style,
      cairowindow::Point position, std::string const& text)
  {
    Text plot_text(convert_to_pixel(position), text);
    add_text(layer, text_style, plot_text);
    return plot_text;
  }

  // Same, but using pixel coordinates.
  [[nodiscard]] Text create_text(boost::intrusive_ptr<Layer> const& layer,
      draw::TextStyle<> const& text_style,
      cairowindow::Pixel position, std::string const& text)
  {
    Text plot_text(position, text);
    add_text(layer, text_style, plot_text);
    return plot_text;
  }

  //--------------------------------------------------------------------------
  // Slider

  [[nodiscard]] Slider create_slider(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Rectangle const& geometry, double start_value, double min_value, double max_value);

  //--------------------------------------------------------------------------

 private:
  void curve_to_lines(boost::intrusive_ptr<Layer> const& layer, Curve const& plot_curve, draw::LineStyle const& line_style);

 public:
  [[nodiscard]] Curve create_curve(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      std::vector<cairowindow::Point>&& points);

  void add_to(boost::intrusive_ptr<Layer> const& layer, bool keep_ratio = false);

  // Called from Window::register_draggable.
  void register_draggable(utils::Badge<Window>, Draggable* draggable,
      std::function<cairowindow::Point (cairowindow::Point const&)>&& restriction)
  {
    ClickableIndex next_index = draggables_.iend();
    draggables_.push_back(draggable);
    draggable->set_index(next_index);
    draggable_restrictions_.emplace_back(std::move(restriction));
  }
  cairowindow::Rectangle update_grabbed(utils::Badge<Window>, ClickableIndex grabbed_point, double pixel_x, double pixel_y);

 private:
  cairowindow::Rectangle axes_geometry(cairowindow::Rectangle const& geometry, double axes_line_width);
  void apply_line_extend(double& x1, double& y1, double& x2, double& y2, LineExtend line_extend);

 public:
#ifdef CWDEBUG
  friend std::ostream& operator<<(std::ostream& os, Plot const* plot_ptr)
  {
    os << "Plot*";
    return os;
  }
#endif
};

} // namespace plot
} // namespace cairowindow
