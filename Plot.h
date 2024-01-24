#pragma once

#include "draw/PlotArea.h"
#include "draw/Text.h"
#include "draw/Point.h"
#include "draw/Circle.h"
#include "draw/Curve.h"
#include "draw/Connector.h"
#include "draw/Arc.h"
#include "draw/Slider.h"
#include "Range.h"
#include "Point.h"
#include "Circle.h"
#include "Text.h"
#include "LinePiece.h"
#include "Line.h"
#include "Curve.h"
#include "Connector.h"
#include "Arc.h"
#include "Draggable.h"
#include "utils/Vector.h"
#include "utils/Badge.h"
#include <boost/intrusive_ptr.hpp>
#include <string>
#include <vector>
#include <functional>
#include <algorithm>
#include "debug.h"
#ifdef CWDEBUG
#include "debug_channel.h"
#endif

namespace cairowindow {
class Window;

namespace plot {
using utils::has_print_on::operator<<;

class Plot;

class Point : public cairowindow::Point, public Draggable
{
 public:
  using cairowindow::Point::Point;
  Point(cairowindow::Point const& point, std::shared_ptr<draw::Point> const& draw_object) :
    cairowindow::Point(point), draw_object_(draw_object) { }

 private:
  friend class Plot;
  friend class cairowindow::Window;
  mutable std::shared_ptr<draw::Point> draw_object_;

  // Implementation of Draggable.
  Rectangle const& geometry() const override { return draw_object_->geometry(); }
  void moved(Plot* plot, cairowindow::Point const& new_position) override;

#ifdef CWDEBUG
  void print_on(std::ostream& os) const override { cairowindow::Point::print_on(os); }
#endif
};

class Circle : public cairowindow::Circle
{
 public:
  using cairowindow::Circle::Circle;
  Circle(cairowindow::Point const& center, double radius, std::shared_ptr<draw::Circle> const& draw_object) :
    cairowindow::Circle(center, radius), draw_object_(draw_object) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Circle> draw_object_;
};

class Text : public cairowindow::Text
{
 public:
  using cairowindow::Text::Text;
  Text(cairowindow::Pixel position, std::string const& text, std::shared_ptr<draw::Text> const& draw_object) :
    cairowindow::Text(position, text), draw_object_(draw_object) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Text> draw_object_;
};

class LinePiece : public cairowindow::LinePiece
{
 public:
  using cairowindow::LinePiece::LinePiece;
  LinePiece(cairowindow::Point const& from, cairowindow::Point const& to, std::shared_ptr<draw::Line> const& draw_object) :
    cairowindow::LinePiece(from, to), draw_object_(draw_object) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

class Line : public cairowindow::Line
{
 public:
  using cairowindow::Line::Line;
  Line(cairowindow::Direction direction, cairowindow::Point point, std::shared_ptr<draw::Line> const& draw_object) :
    cairowindow::Line(direction, point), draw_object_(draw_object) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Line> draw_object_;
};

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

class Connector : public cairowindow::Connector
{
 public:
  using cairowindow::Connector::Connector;
  Connector(cairowindow::Point const& from, cairowindow::Point const& to,
      ArrowHeadShape head_from, ArrowHeadShape head_to,
      std::shared_ptr<draw::Connector> const& draw_object) :
    cairowindow::Connector(from, to, head_from, head_to), draw_object_(draw_object) { }

  Connector(cairowindow::Point const& from, cairowindow::Point const& to,
      ArrowHeadShape head_to, std::shared_ptr<draw::Connector> const& draw_object) :
    cairowindow::Connector(from, to, head_to), draw_object_(draw_object) { }

  Connector(cairowindow::Point const& from, cairowindow::Point const& to,
      std::shared_ptr<draw::Connector> const& draw_object) :
    cairowindow::Connector(from, to), draw_object_(draw_object) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Connector> draw_object_;
};

class Arc : public cairowindow::Arc
{
 public:
  using cairowindow::Arc::Arc;
  Arc(Point const& center, double start_angle, double end_angle, double radius,
      std::shared_ptr<draw::Arc> const& draw_object) :
    cairowindow::Arc(center, start_angle, end_angle, radius), draw_object_(draw_object) { }

  Arc(Point const& center, Direction const& start, Direction const& end, double radius,
      std::shared_ptr<draw::Arc> const& draw_object) :
    cairowindow::Arc(center, start, end, radius), draw_object_(draw_object) { }

  Arc(Line const& line1, Line const& line2, double radius,
      std::shared_ptr<draw::Arc> const& draw_object) :
    cairowindow::Arc(line1, line2, radius), draw_object_(draw_object) { }

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Arc> draw_object_;
};

class Slider
{
 private:
  Rectangle geometry_;
  double min_value_;      // The value corresponding to position 0.
  double max_value_;      // The value corresponding to position 1.

 public:
  Slider(Rectangle geometry, double min_value, double max_value) :
    geometry_(geometry), min_value_(min_value), max_value_(max_value) { }

  Rectangle const& geometry() const { return geometry_; }
  double value() const { return min_value_ + draw_object_->rel_value() * (max_value_ - min_value_); }
  double min_value() const { return min_value_; }
  double max_value() const { return max_value_; }

  void set_value(double value);

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Slider> draw_object_;
};

enum class LineExtend
{
  none = 0,
  from = 1,
  to = 2,
  both = from|to
};

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
  Plot(Rectangle const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, TitleStyle title_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width), plot_area_style),
    title_(std::make_shared<draw::Text>(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_style.offset, title_style)) { }

  Plot(Rectangle const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, TitleStyle title_style,
      std::string xlabel, XLabelStyle xlabel_style, std::string ylabel, YLabelStyle ylabel_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width), plot_area_style),
    title_(std::make_shared<draw::Text>(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_style.offset, title_style)),
    xlabel_(std::make_shared<draw::Text>(xlabel, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() + plot_area_.geometry().height() + XLabelStyleDefaults::offset, xlabel_style)),
    ylabel_(std::make_shared<draw::Text>(ylabel, plot_area_.geometry().offset_x() - YLabelStyleDefaults::offset,
        plot_area_.geometry().offset_y() + 0.5 * plot_area_.geometry().height(), ylabel_style)) { }

  Plot(Rectangle const& geometry, draw::PlotAreaStyle plot_area_style) :
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

  // Create and draw a point on layer at x,y using point_style.
  [[nodiscard]] Point create_point(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point const& point, draw::PointStyle const& point_style);

  // Create and draw a circle on layer with center and radius using circle_style.
  [[nodiscard]] Circle create_circle(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point const& center, double radius, draw::CircleStyle const& circle_style);

  // Same as above but use line_style (no fill_color).
  [[nodiscard]] Circle create_circle(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point const& center, double radius, draw::LineStyle const& line_style)
  {
    return create_circle(layer, center, radius, draw::CircleStyle{.line_color = line_style.line_color, .line_width = line_style.line_width});
  }

  // Create and drw an arc on layer width center, radius and start- and end_angle, using arc_style.
  [[nodiscard]] Arc create_arc(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point const& center, double start_angle, double end_angle, double radius, draw::ArcStyle const& arc_style);

  // Create and draw text on layer at position using text_style.
  [[nodiscard]] Text create_text(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point position, std::string const& text, draw::TextStyle<> const& text_style);

  // Same, but using pixel coordinates.
  [[nodiscard]] Text create_text(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Pixel position, std::string const& text, draw::TextStyle<> const& text_style);

  // Create and draw a line piece between points from and to using line_style and line_extend.
  [[nodiscard]] LinePiece create_line(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point const& from, cairowindow::Point const& to, draw::LineStyle const& line_style,
      LineExtend line_extend = LineExtend::none);

  // Create and draw a line piece between points from and to using line_style and line_extend.
  [[nodiscard]] LinePiece create_line(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::LinePiece const& line_piece, draw::LineStyle const& line_style,
      LineExtend line_extend = LineExtend::none)
  {
    return create_line(layer, line_piece.from(), line_piece.to(), line_style, line_extend);
  }

  // Create and draw a line through point in direction using line_style.
  [[nodiscard]] Line create_line(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point const& point, cairowindow::Direction const& direction, draw::LineStyle const& line_style);

  // Create and draw a line according to Line using line_style.
  [[nodiscard]] Line create_line(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Line const& line, draw::LineStyle const& line_style)
  {
    return create_line(layer, line.point(), line.direction(), line_style);
  }

  // Create and draw a connector from point to point using line_style and fill_color for the arrow heads if appropriate.
  [[nodiscard]] Connector create_connector(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point const& from, cairowindow::Point const& to,
      Connector::ArrowHeadShape arrow_head_shape_from, Connector::ArrowHeadShape arrow_head_shape_to,
      draw::LineStyle const& line_style, Color fill_color = color::white);

  [[nodiscard]] Connector create_connector(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point const& from, cairowindow::Point const& to,
      Connector::ArrowHeadShape arrow_head_shape_to,
      draw::LineStyle const& line_style, Color fill_color = color::white)
  {
    return create_connector(layer, from, to, Connector::no_arrow, arrow_head_shape_to, line_style, fill_color);
  }

  [[nodiscard]] Connector create_connector(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Point const& from, cairowindow::Point const& to,
      draw::LineStyle const& line_style, Color fill_color = color::white)
  {
    return create_connector(layer, from, to, Connector::no_arrow, Connector::open_arrow, line_style, fill_color);
  }

  [[nodiscard]] Slider create_slider(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Rectangle const& geometry, double start_value, double min_value, double max_value);

  [[nodiscard]] Curve create_curve(boost::intrusive_ptr<Layer> const& layer,
      std::vector<cairowindow::Point>&& points, draw::LineStyle const& line_style);

  // Same, but pass an already existing plot_point. If plot_point was added before
  // then it will be removed from the plot and then added to the new coordinates.
  void add_point(boost::intrusive_ptr<Layer> const& layer,
      Point const& plot_point, draw::PointStyle const& point_style);

  void add_text(boost::intrusive_ptr<Layer> const& layer,
      Text const& text, draw::TextStyle<> const& text_style);

  void add_line(boost::intrusive_ptr<Layer> const& layer, LinePiece const& plot_line_piece,
      draw::LineStyle const& line_style, LineExtend line_extend = LineExtend::none);

  void add_connector(boost::intrusive_ptr<Layer> const& layer, Connector const& plot_connector,
      draw::LineStyle const& line_style, Color fill_color = color::white);

  void add_arc(boost::intrusive_ptr<Layer> const& layer, Arc const& plot_arc,
      draw::ArcStyle const& arc_style);

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
  Rectangle update_grabbed(utils::Badge<Window>, ClickableIndex grabbed_point, double pixel_x, double pixel_y);

 private:
  Rectangle axes_geometry(Rectangle const& geometry, double axes_line_width);
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
