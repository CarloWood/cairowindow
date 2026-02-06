#pragma once

#include "CoordinateMapper.h"
#include "Range.h"
#include "Circle.h"
#include "Text.h"
#include "Curve.h"
#include "Arc.h"
#include "BezierCurve.h"
#include "BezierFitter.h"
#include "Draggable.h"
#include "plot/Connector.h"
#include "plot/Point.h"
#include "plot/Rectangle.h"
#include "plot/LinePiece.h"
#include "plot/Line.h"
#include "cs/Point.h"
#include "draw/PlotArea.h"
#include "draw/Text.h"
#include "draw/Point.h"
#include "draw/Circle.h"
#include "draw/Curve.h"
#include "draw/Connector.h"
#include "draw/Rectangle.h"
#include "draw/Arc.h"
#include "draw/BezierCurve.h"
#include "draw/BezierFitter.h"
#include "draw/Slider.h"
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
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Plot;

//--------------------------------------------------------------------------
// Slider

class Slider
{
 private:
  cairowindow::Geometry geometry_;
  double min_value_;      // The value corresponding to position 0.
  double max_value_;      // The value corresponding to position 1.

 public:
  Slider(cairowindow::Geometry geometry, double min_value, double max_value) :
    geometry_(geometry), min_value_(min_value), max_value_(max_value) { }

  cairowindow::Geometry const& geometry() const { return geometry_; }
  double value() const { return min_value_ + draw_object_->rel_value() * (max_value_ - min_value_); }
  double min_value() const { return min_value_; }
  double max_value() const { return max_value_; }

  void set_value(double value);

 public:
  mutable std::shared_ptr<draw::Slider> draw_object_;
};

//--------------------------------------------------------------------------

class Curve : public cairowindow::Curve
{
 public:
  Curve(std::vector<cairowindow::Point> const& points, std::shared_ptr<draw::Curve> const& draw_object) :
    cairowindow::Curve(points), draw_object_(draw_object) { }
  Curve(std::vector<cairowindow::Point>&& points, std::shared_ptr<draw::Curve> const& draw_object) :
    cairowindow::Curve(std::move(points)), draw_object_(draw_object) { }

  explicit Curve(cairowindow::Curve const& curve) : cairowindow::Curve(curve) { }
  using cairowindow::Curve::Curve;

 public:
  mutable std::shared_ptr<draw::Curve> draw_object_;
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

struct TitleStyleDefaults : draw::TextStyleParamsDefault
{
  static constexpr draw::TextPosition position = draw::centered;
  static constexpr double font_size = 24.0;
};

struct LabelStyleDefaults : draw::TextStyleParamsDefault
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

} // namespace plot
namespace draw {

#define cairowindow_PlotTitle_FOREACH_MEMBER cairowindow_TextBase_FOREACH_MEMBER
#define cairowindow_PlotTitle_FOREACH_STYLE_MEMBER cairowindow_TextBase_FOREACH_MEMBER
#define cairowindow_XLabel_FOREACH_MEMBER cairowindow_TextBase_FOREACH_MEMBER
#define cairowindow_XLabel_FOREACH_STYLE_MEMBER cairowindow_TextBase_FOREACH_MEMBER
#define cairowindow_YLabel_FOREACH_MEMBER cairowindow_TextBase_FOREACH_MEMBER
#define cairowindow_YLabel_FOREACH_STYLE_MEMBER cairowindow_TextBase_FOREACH_MEMBER
#define cairowindow_Label_FOREACH_MEMBER cairowindow_TextBase_FOREACH_MEMBER
#define cairowindow_Label_FOREACH_STYLE_MEMBER cairowindow_TextBase_FOREACH_MEMBER

DECLARE_STYLE_WITH_BASE(PlotTitle, TextBase, plot::TitleStyleDefaults)
DECLARE_STYLE_WITH_BASE(XLabel, TextBase, plot::XLabelStyleDefaults)
DECLARE_STYLE_WITH_BASE(YLabel, TextBase, plot::YLabelStyleDefaults)
DECLARE_STYLE_WITH_BASE(Label, TextBase, plot::LabelStyleDefaults)

#undef cairowindow_PlotTitle_FOREACH_MEMBER
#undef cairowindow_PlotTitle_FOREACH_STYLE_MEMBER
#undef cairowindow_XLabel_FOREACH_MEMBER
#undef cairowindow_XLabel_FOREACH_STYLE_MEMBER
#undef cairowindow_YLabel_FOREACH_MEMBER
#undef cairowindow_YLabel_FOREACH_STYLE_MEMBER
#undef cairowindow_Label_FOREACH_MEMBER
#undef cairowindow_Label_FOREACH_STYLE_MEMBER

} // namespace draw
namespace plot {

class Plot : public CoordinateMapper<csid::plot>
{
 public:
  using Vector = math::Vector<2>;

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

 public:
  Plot(cairowindow::Geometry const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, draw::PlotTitleStyle title_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width()), plot_area_style),
    title_(std::make_shared<draw::Text>(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_style.offset(), title_style)) { }

  Plot(cairowindow::Geometry const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, draw::PlotTitleStyle title_style,
      std::string xlabel, draw::XLabelStyle xlabel_style, std::string ylabel, draw::YLabelStyle ylabel_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width()), plot_area_style),
    title_(std::make_shared<draw::Text>(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_style.offset(), title_style)),
    xlabel_(std::make_shared<draw::Text>(xlabel, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() + plot_area_.geometry().height() + XLabelStyleDefaults::offset, xlabel_style)),
    ylabel_(std::make_shared<draw::Text>(ylabel, plot_area_.geometry().offset_x() - YLabelStyleDefaults::offset,
        plot_area_.geometry().offset_y() + 0.5 * plot_area_.geometry().height(), ylabel_style)) { }

  Plot(cairowindow::Geometry const& geometry, draw::PlotAreaStyle plot_area_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width()), plot_area_style) { }

  void set_range(int axis, Range range)
  {
    DoutEntering(dc::cairowindow, "Plot::set_range(" << axis << ", " << range << ") [" << this << "]");
    range_[axis] = range;
    range_ticks_[axis] = draw::PlotArea::calculate_range_ticks(range_[axis]);
    Dout(dc::cairowindow, "range_[" << axis << "] = " << range_[axis] << "; range_ticks_[" << axis << "] = " << range_ticks_[axis]);
    update_plot_transform_pixels();
  }

  void set_xrange(Range x_range) { set_range(x_axis, x_range); }
  void set_yrange(Range y_range) { set_range(y_axis, y_range); }

  cairowindow::Point clamp_to_plot_area(cairowindow::Point const& point) const
  {
    return {std::clamp(point.x(), range_[x_axis].min(), range_[x_axis].max()),
            std::clamp(point.y(), range_[y_axis].min(), range_[y_axis].max())};
  }

  Range const& xrange() const { return range_[x_axis]; }
  Range const& yrange() const { return range_[y_axis]; }

  cairowindow::Rectangle viewport() const
  {
    return {range_[x_axis].min(), range_[y_axis].min(), range_[x_axis].size(), range_[y_axis].size()};
  }

  double convert_x(double x) const;
  double convert_y(double y) const;
  Pixel convert_to_pixel(cairowindow::Point const& point) const;

  double convert_from_pixel_x(double pixel_x) const;
  double convert_from_pixel_y(double pixel_y) const;
  cairowindow::Point convert_from_pixel(Pixel const& pixel) const;

  double convert_horizontal_offset_from_pixel(double pixel_offset_x) const;
  double convert_vertical_offset_from_pixel(double pixel_offset_y) const;

  void convert_to_pixels(cairowindow::Point const* data_in, Pixel* data_out, std::size_t size);

  template<std::size_t size>
  [[gnu::always_inline]] void convert_to_pixels(std::array<cairowindow::Point, size> const& in, std::array<Pixel, size>& out)
  {
    convert_to_pixels(in.data(), out.data(), size);
  }

  //--------------------------------------------------------------------------
  // Line

  // Add and draw plot_line using line_style.
  void add_line(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      Line const& plot_line);

 public:
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
  // LinePiece

  // Add and draw plot_line_piece on layer using line_style and line_extend.
  void add_line_piece(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style, LineExtend line_extend,
      LinePiece const& plot_line_piece);

 public:
  // Create and draw a line piece between points from and to using line_style and line_extend.
  template<typename... Args>
  [[nodiscard]] LinePiece create_line_piece(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style, LineExtend line_extend,
      Args&&... args)
    requires requires(Args&&... args) { LinePiece{std::forward<Args>(args)...}; }
  {
    LinePiece plot_line_piece(std::forward<Args>(args)...);
    add_line_piece(layer, line_style, line_extend, plot_line_piece);
    return plot_line_piece;
  }

  // Same, but without a line_extend.
  template<typename... Args>
  [[nodiscard]] LinePiece create_line_piece(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      Args&&... args)
    requires requires(Args&&... args) { LinePiece{std::forward<Args>(args)...}; }
  {
    LinePiece plot_line_piece(std::forward<Args>(args)...);
    add_line_piece(layer, line_style, LineExtend::none, plot_line_piece);
    return plot_line_piece;
  }

  //--------------------------------------------------------------------------
  // Connector

 private:
  template<typename... Args>
  [[nodiscard]] Connector create_connector_helper(
      Connector::ArrowHeadShape& default_arrow_head_shape_from,
      Connector::ArrowHeadShape& default_arrow_head_shape_to,
      boost::intrusive_ptr<Layer> const& layer,
      draw::ConnectorStyle const& connector_style, std::tuple<Args...>&& args)
  {
    // Strip Connector::ArrowHeadShape arguments from args... and apply them to the defaults.
    if constexpr (std::is_same_v<std::tuple_element_t<0, std::tuple<Args...>>, typename Connector::ArrowHeadShape>)
    {
      if constexpr (std::is_same_v<std::tuple_element_t<1, std::tuple<Args...>>, typename Connector::ArrowHeadShape>)
        default_arrow_head_shape_from = std::get<0>(args);
      else
        default_arrow_head_shape_to = std::get<0>(args);
      return create_connector_helper(default_arrow_head_shape_from, default_arrow_head_shape_to,
          layer, connector_style, tuple_tail(std::move(args)));
    }
    else
    {
      // The defaults are now set. Construct a Connector from the remaining arguments and the default arrow head shapes.
      Connector plot_connector = std::apply([&](auto&&... unpacked_args) -> Connector {
        return {std::forward<decltype(unpacked_args)>(unpacked_args)..., default_arrow_head_shape_from, default_arrow_head_shape_to};
      }, std::move(args));

      add_connector(layer, connector_style, plot_connector);
      return plot_connector;
    }
  }

 public:
  // Create and draw a connector on layer, using args... and line_style.
  // Args can optionally start with zero, one or two ArrowHeadShape arguments.
  template<typename... Args>
  [[nodiscard]] Connector create_connector(boost::intrusive_ptr<Layer> const& layer,
      draw::ConnectorStyle const& connector_style, Args&&... args)
  {
    Connector::ArrowHeadShape arrow_head_shape_from = Connector::no_arrow;
    Connector::ArrowHeadShape arrow_head_shape_to = Connector::open_arrow;
    return create_connector_helper(arrow_head_shape_from, arrow_head_shape_to, layer,
        connector_style, std::make_tuple(std::forward<Args>(args)...));
  }

  //--------------------------------------------------------------------------
  // Rectangle

 public:
  // Create and draw a rectangle on layer, using args... and rectangle_style.
  template<typename... Args>
  [[nodiscard]] Rectangle create_rectangle(boost::intrusive_ptr<Layer> const& layer,
      draw::RectangleStyle const& rectangle_style, Args&&... args)
  {
    Rectangle plot_rectangle(std::forward<Args>(args)...);
    add_rectangle(layer, rectangle_style, plot_rectangle);
    return plot_rectangle;
  }

  //--------------------------------------------------------------------------
  // Circle

 public:
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

 public:
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

 private:
  void add_bezier_curve( // Do not pass a cairowindow::BezierCurve to this function! It must be a plot::BezierCurve.
      boost::intrusive_ptr<Layer> const& layer,
      draw::BezierCurveStyle const& bezier_style,
      BezierCurve&& plot_bezier_curve);

 public:
  // Add and draw plot_bezier_curve on layer using bezier_style.
  // The difference with the previous one is that in this case it is assumed that the x and y coordinates
  // of the Bezier curve are in pixels already, and no conversion will take place.
  void add_bezier_curve_in_px(boost::intrusive_ptr<Layer> const& layer,
      draw::BezierCurveStyle const& bezier_style,
      BezierCurve const& plot_bezier_curve_in_px);

 private:
  void add_bezier_curve_in_px( // Do not pass a cairowindow::BezierCurve to this function! It must be a plot::BezierCurve.
      boost::intrusive_ptr<Layer> const& layer,
      draw::BezierCurveStyle const& bezier_style,
      BezierCurve&& plot_bezier_curve_in_px);

 public:
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

  // Create and draw a Bezier curve on layer using bezier_style.
  // The difference with the previous one is that it is assumed the Args passed are in pixels:
  // the constructed BezierCurve is added using add_bezier_curve_in_px.
  template<typename... Args>
  [[nodiscard]] BezierCurve create_bezier_curve_in_px(boost::intrusive_ptr<Layer> const& layer,
      draw::BezierCurveStyle const& bezier_style,
      Args&&... args)
  {
    BezierCurve plot_bezier_curve_in_px(std::forward<Args>(args)...);
    add_bezier_curve_in_px(layer, bezier_style, plot_bezier_curve_in_px);
    return plot_bezier_curve_in_px;
  }

  //--------------------------------------------------------------------------
  // Curve

  void add_bezier_fitter(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
      draw::PointStyle const& point_style,
#endif
      BezierFitter const& plot_bezier_fitter);

 private:
  void add_bezier_fitter( // Do not pass a cairowindow::BezierFitter to this function! It must be a plot::BezierFitter.
      boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
      draw::PointStyle const& point_style,
#endif
      BezierFitter&& plot_bezier_fitter);

 public:
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
  void add_bezier_fitter(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      BezierFitter const& plot_bezier_fitter)
  {
    add_bezier_fitter(layer, line_style, draw::PointStyle{}, plot_bezier_fitter);
  }

 private:
  void add_bezier_fitter( // Do not pass a cairowindow::BezierFitter to this function! It must be a plot::BezierFitter.
      boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      BezierFitter&& plot_bezier_fitter);
#endif

  [[nodiscard]] BezierFitter create_bezier_fitter(boost::intrusive_ptr<Layer> const& layer,
      draw::LineStyle const& line_style,
      cairowindow::BezierFitter&& bezier_fitter)
  {
    BezierFitter plot_bezier_fitter(std::move(bezier_fitter));
    add_bezier_fitter(layer, line_style, plot_bezier_fitter);
    return plot_bezier_fitter;
  }

  //--------------------------------------------------------------------------
  // Text

  // Add and draw plot_text on layer using text_style.
  void add_text(boost::intrusive_ptr<Layer> const& layer,
      draw::TextStyle const& text_style,
      Text const& plot_text);

 private:
  void add_text( // Do not pass a cairowindow::Text to this function! It must be a plot::Text.
      boost::intrusive_ptr<Layer> const& layer,
      draw::TextStyle const& text_style,
      Text&& plot_text);

 public:
  // Create and draw text on layer at position using text_style.
  [[nodiscard]] Text create_text(boost::intrusive_ptr<Layer> const& layer,
      draw::TextStyle const& text_style,
      cairowindow::Point position, std::string const& text)
  {
    Text plot_text(convert_to_pixel(position), text);
    add_text(layer, text_style, plot_text);
    return plot_text;
  }

  // Same, but using pixel coordinates.
  [[nodiscard]] Text create_text(boost::intrusive_ptr<Layer> const& layer,
      draw::TextStyle const& text_style,
      cairowindow::Pixel position, std::string const& text)
  {
    Text plot_text(position, text);
    add_text(layer, text_style, plot_text);
    return plot_text;
  }

  //--------------------------------------------------------------------------
  // Slider

  [[nodiscard]] Slider create_slider(boost::intrusive_ptr<Layer> const& layer,
      cairowindow::Geometry const& geometry, double start_value, double min_value, double max_value);

  //--------------------------------------------------------------------------

 private:
  void update_plot_transform_pixels();
  void curve_to_bezier_curves(boost::intrusive_ptr<Layer> const& layer,
      Curve const& plot_curve, draw::BezierCurveStyle const& bezier_curve_style);

 public:
  [[nodiscard]] Curve create_curve(boost::intrusive_ptr<Layer> const& layer,
      draw::BezierCurveStyle const& line_style,
      std::vector<cairowindow::Point>&& points);

  cairowindow::Geometry const& geometry() const override { return plot_area_.geometry(); }
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
  // Called from Window::update_grabbed and Window::move_draggable.
  cairowindow::Geometry update_draggable(utils::Badge<Window>, ClickableIndex draggable_index, cairowindow::Point const& new_position);
  cairowindow::Geometry update_grabbed(utils::Badge<Window>, ClickableIndex grabbed_point, double pixel_x, double pixel_y);
  void apply_restrictions(utils::Badge<Window>, ClickableIndex clickable_index, cairowindow::Point& new_position);

 private:
  cairowindow::Geometry axes_geometry(cairowindow::Geometry const& geometry, double axes_line_width);

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
