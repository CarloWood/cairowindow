#pragma once

#include "Printable.h"
#include "Text.h"
#include "plot/Arc.h"
#include "plot/Circle.h"
#include "plot/Connector.h"
#include "plot/Rectangle.h"
#include "plot/LinePiece.h"
#include "plot/Line.h"
#include "plot/Point.h"
#include "draw/Arc.h"
#include "draw/Circle.h"
#include "draw/Connector.h"
#include "draw/Rectangle.h"
#include "draw/Polyline.h"
#include "draw/Line.h"
#include "draw/Point.h"
#include "draw/Text.h"
#include "math/cs/Point.h"
#include "math/cs/TransformOperators.h"
#include "math/Hyperblock.h"
#include "math/Transform.h"
#include <cmath>
#include <functional>
#include <tuple>
#include <type_traits>
#include <utility>

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

namespace cairowindow {

// Type passed to CoordinateMapper<cs>::add_clipped_line_piece; determines if a line piece is extended beyond the given end points.
enum class LineExtend
{
  none = 0,
  from = 1,
  to = 2,
  both = from | to
};

//-----------------------------------------------------------------------------
// CoordinateMapper
//
// A small utility base class that maps logical coordinate system `cs` to pixels
// using a math::Transform<cs, csid::pixels>. It offers shared implementation for
// drawing and dragging plot::cs::Point<cs> instances.
//
template<CS cs>
class CoordinateMapper : public Printable
{
 protected:
  using LayerPtr = boost::intrusive_ptr<Layer>;
  using PointHandle = plot::cs::Point<cs>;
  using RectangleHandle = plot::cs::Rectangle<cs>;
  using LineHandle = plot::cs::Line<cs>;
  using LinePieceHandle = plot::cs::LinePiece<cs>;
  using CircleHandle = plot::cs::Circle<cs>;
  using ArcHandle = plot::cs::Arc<cs>;
  using ConnectorHandle = plot::cs::Connector<cs>;
  using TextHandle = plot::Text;

 protected:
  math::Transform<cs, csid::pixels> cs_transform_pixels_;

 public:
  CoordinateMapper() = default;
  explicit CoordinateMapper(math::Transform<cs, csid::pixels> const& cs_transform_pixels) :
    cs_transform_pixels_(cs_transform_pixels) { }

  math::Transform<cs, csid::pixels> const& cs_transform_pixels() const { return cs_transform_pixels_; }
  void set_cs_transform_pixels(math::Transform<cs, csid::pixels> const& cs_transform_pixels)
  {
    cs_transform_pixels_ = cs_transform_pixels;
  }

  //--------------------------------------------------------------------------
  // Point

  // Add and draw plot_point_cs on layer using point_style.
  void add_point(LayerPtr const& layer, draw::PointStyle const& point_style, PointHandle const& plot_point_cs);

  // Create and draw a point on layer using point_style.
  [[nodiscard]] PointHandle create_point(LayerPtr const& layer, draw::PointStyle const& point_style, cs::Point<cs> const& point_cs)
  {
    PointHandle plot_point_cs{point_cs};
    add_point(layer, point_style, plot_point_cs);
    return plot_point_cs;
  }

  // Called by Window::update_grabbed through the lambda defined in Window::register_draggable<cs> when a Draggable plot_point_cs was moved to (pixel_x, pixel_y).
  Geometry update_grabbed(plot::cs::Point<cs>* plot_point_cs, double pixel_x, double pixel_y,
      std::function<cs::Point<cs> (cs::Point<cs> const&)> const& restriction);

  //--------------------------------------------------------------------------
  // Text

  // Add and draw plot_text on layer using text_style.
  void add_text(LayerPtr const& layer, draw::TextStyle const& text_style, TextHandle const& plot_text);

  // Create and draw text on layer at position using text_style.
  [[nodiscard]] TextHandle create_text(LayerPtr const& layer, draw::TextStyle const& text_style,
      cs::Point<cs> const& position, std::string const& text)
  {
    cs::Point<csid::pixels> const position_pixels = position * cs_transform_pixels_;
    TextHandle plot_text(position_pixels, text);
    add_text(layer, text_style, plot_text);
    return plot_text;
  }

  // Same, but using pixel coordinates.
  [[nodiscard]] TextHandle create_text(LayerPtr const& layer, draw::TextStyle const& text_style,
      Pixel position, std::string const& text)
  {
    TextHandle plot_text(position, text);
    add_text(layer, text_style, plot_text);
    return plot_text;
  }

  //--------------------------------------------------------------------------
  // Rectangle

  // Add and draw plot_rectangle_cs on layer using rectangle_style. Handles rotation by drawing a polyline.
  void add_rectangle(LayerPtr const& layer, draw::RectangleStyle const& rectangle_style, RectangleHandle const& plot_rectangle_cs);

  // Create and draw a rectangle on layer, using args... and rectangle_style.
  template<typename... Args>
  [[nodiscard]] RectangleHandle create_rectangle(LayerPtr const& layer, draw::RectangleStyle const& rectangle_style, Args&&... args)
    requires requires(Args&&... args) { RectangleHandle{std::forward<Args>(args)...}; }
  {
    RectangleHandle plot_rectangle_cs(std::forward<Args>(args)...);
    add_rectangle(layer, rectangle_style, plot_rectangle_cs);
    return plot_rectangle_cs;
  }

  //--------------------------------------------------------------------------
  // Circle

  // Add and draw plot_circle_cs, with its center and radius in cs coordinates, on layer using circle_style.
  void add_circle(LayerPtr const& layer, draw::CircleStyle const& circle_style, CircleHandle const& plot_circle_cs);

  // Create and draw a circle on layer with center and radius using circle_style.
  template<typename... Args>
  [[nodiscard]] CircleHandle create_circle(LayerPtr const& layer, draw::CircleStyle const& circle_style, Args&&... args)
    requires requires(Args&&... args) { CircleHandle{std::forward<Args>(args)...}; }
  {
    CircleHandle plot_circle_cs(std::forward<Args>(args)...);
    add_circle(layer, circle_style, plot_circle_cs);
    return plot_circle_cs;
  }

  // Same as above but use line_style (no fill_color).
  template<typename... Args>
  [[nodiscard]] CircleHandle create_circle(LayerPtr const& layer, draw::LineStyle const& line_style, Args&&... args)
  {
    return create_circle(layer, draw::CircleStyle({.line_color = line_style.line_color(), .line_width = line_style.line_width()}),
        std::forward<Args>(args)...);
  }

  //--------------------------------------------------------------------------
  // Arc

  // Add and draw plot_arc_cs on layer using arc_style.
  void add_arc(LayerPtr const& layer, draw::ArcStyle const& arc_style, ArcHandle const& plot_arc_cs);

  // Create and draw an arc on layer with center, radius and start- and end_angle, using arc_style.
  template<typename... Args>
  [[nodiscard]] ArcHandle create_arc(LayerPtr const& layer, draw::ArcStyle const& arc_style, Args&&... args)
    requires requires(Args&&... args) { ArcHandle{std::forward<Args>(args)...}; }
  {
    ArcHandle plot_arc_cs(std::forward<Args>(args)...);
    add_arc(layer, arc_style, plot_arc_cs);
    return plot_arc_cs;
  }

  //--------------------------------------------------------------------------
  // Line (infinite, clipped)

  // Add and draw plot_line_cs on layer using line_style, clipped to clip_rectangle_cs.
  void add_clipped_line(LayerPtr const& layer, draw::LineStyle const& line_style,
      LineHandle const& plot_line_cs, math::Hyperblock<2> const& clip_rectangle_cs);

  //--------------------------------------------------------------------------
  // LinePiece

  // Add and draw plot_line_piece_cs on layer using line_style and line_extend, clipped to clip_rectangle_cs.
  void add_clipped_line_piece(LayerPtr const& layer, draw::LineStyle const& line_style, LineExtend line_extend,
      LinePieceHandle const& plot_line_piece_cs, math::Hyperblock<2> const& clip_rectangle_cs);

  //--------------------------------------------------------------------------
  // Connector

  // Add and draw plot_connector_cs on layer using connector_style.
  void add_connector(LayerPtr const& layer, draw::ConnectorStyle const& connector_style, ConnectorHandle const& plot_connector_cs);

 private:
  template<typename... Args>
  [[nodiscard]] ConnectorHandle create_connector_helper(
      typename ConnectorHandle::ArrowHeadShape& default_arrow_head_shape_from,
      typename ConnectorHandle::ArrowHeadShape& default_arrow_head_shape_to,
      LayerPtr const& layer,
      draw::ConnectorStyle const& connector_style, std::tuple<Args...>&& args)
  {
    // Strip ConnectorHandle::ArrowHeadShape arguments from args... and apply them to the defaults.
    if constexpr (std::is_same_v<std::tuple_element_t<0, std::tuple<Args...>>, typename ConnectorHandle::ArrowHeadShape>)
    {
      if constexpr (std::is_same_v<std::tuple_element_t<1, std::tuple<Args...>>, typename ConnectorHandle::ArrowHeadShape>)
        default_arrow_head_shape_from = std::get<0>(args);
      else
        default_arrow_head_shape_to = std::get<0>(args);
      return create_connector_helper(default_arrow_head_shape_from, default_arrow_head_shape_to,
          layer, connector_style, tuple_tail(std::move(args)));
    }
    else
    {
      // The defaults are now set. Construct a ConnectorHandle from the remaining arguments and the default arrow head shapes.
      ConnectorHandle plot_connector_cs = std::apply([&](auto&&... unpacked_args) -> ConnectorHandle {
        return {std::forward<decltype(unpacked_args)>(unpacked_args)..., default_arrow_head_shape_from, default_arrow_head_shape_to};
      }, std::move(args));

      add_connector(layer, connector_style, plot_connector_cs);
      return plot_connector_cs;
    }
  }

 public:
  // Create and draw a connector on layer, using args... and connector_style.
  // Args can optionally start with zero, one or two ArrowHeadShape arguments.
  template<typename... Args>
  [[nodiscard]] ConnectorHandle create_connector(LayerPtr const& layer, draw::ConnectorStyle const& connector_style, Args&&... args)
  {
    typename ConnectorHandle::ArrowHeadShape arrow_head_shape_from = ConnectorHandle::no_arrow;
    typename ConnectorHandle::ArrowHeadShape arrow_head_shape_to = ConnectorHandle::open_arrow;
    return create_connector_helper(arrow_head_shape_from, arrow_head_shape_to, layer,
        connector_style, std::make_tuple(std::forward<Args>(args)...));
  }
};

//----------------------------------------------------------------------------
// PointHandle

template<CS cs>
void CoordinateMapper<cs>::add_point(LayerPtr const& layer, draw::PointStyle const& point_style, PointHandle const& plot_point_cs)
{
  // Convert from the coordinate-system space to pixels using the transform.
  cs::Point<csid::pixels> point_pixels = plot_point_cs * cs_transform_pixels_;

  plot_point_cs.create_draw_object({}, point_pixels.x(), point_pixels.y(), point_style);
  draw_layer_region_on(layer, plot_point_cs.draw_object());
}

template<CS cs>
Geometry CoordinateMapper<cs>::update_grabbed(plot::cs::Point<cs>* plot_point_cs, double pixel_x, double pixel_y,
    std::function<cs::Point<cs> (cs::Point<cs> const&)> const& restriction)
{
  cs::Point<csid::pixels> const new_position_pixels{pixel_x, pixel_y};
  cs::Point<cs> new_position_cs = new_position_pixels * cs_transform_pixels_.inverse();

  if (restriction)
    new_position_cs = restriction(new_position_cs);

  *plot_point_cs = new_position_cs;

  auto const& draw_object = plot_point_cs->draw_object();
  add_point(draw_object->layer(), draw_object->point_style(), *plot_point_cs);

  return plot_point_cs->geometry();
}

//----------------------------------------------------------------------------
// TextHandle

template<CS cs>
void CoordinateMapper<cs>::add_text(LayerPtr const& layer, draw::TextStyle const& text_style, TextHandle const& plot_text)
{
  cairowindow::Pixel const position = plot_text.position();
  std::string const& text = plot_text.text();

  plot_text.template create_draw_object<cs>({}, text, position.x(), position.y(), text_style);
  draw_layer_region_on(layer, plot_text.draw_object());
}

//----------------------------------------------------------------------------
// RectangleHandle

template<CS cs>
void CoordinateMapper<cs>::add_rectangle(LayerPtr const& layer, draw::RectangleStyle const& rectangle_style, RectangleHandle const& plot_rectangle_cs)
{
  // Use a fast axis-aligned pixel rectangle when the transform maps the basis vectors to axis-parallel directions.
  // Otherwise, draw a rotated/sheared rectangle outline as a polyline.
  auto const [vx_x, vx_y] = cs_transform_pixels_.map_vector(1.0, 0.0);
  auto const [vy_x, vy_y] = cs_transform_pixels_.map_vector(0.0, 1.0);

  double constexpr eps = 1e-12;
  bool const vx_parallel_x = std::abs(vx_y) < eps;
  bool const vx_parallel_y = std::abs(vx_x) < eps;
  bool const vy_parallel_x = std::abs(vy_y) < eps;
  bool const vy_parallel_y = std::abs(vy_x) < eps;
  bool const basis_axis_aligned = (vx_parallel_x || vx_parallel_y) && (vy_parallel_x || vy_parallel_y) &&
    ((vx_parallel_x && vy_parallel_y) || (vx_parallel_y && vy_parallel_x));

  if (basis_axis_aligned)
  {
    cs::Point<cs> const p1_cs(plot_rectangle_cs.offset_x(), plot_rectangle_cs.offset_y());
    cs::Point<cs> const p2_cs(plot_rectangle_cs.offset_x() + plot_rectangle_cs.width(), plot_rectangle_cs.offset_y() + plot_rectangle_cs.height());

    cs::Point<csid::pixels> const p1_pixels = p1_cs * cs_transform_pixels_;
    cs::Point<csid::pixels> const p2_pixels = p2_cs * cs_transform_pixels_;

    plot_rectangle_cs.create_draw_object({}, p1_pixels.x(), p1_pixels.y(), p2_pixels.x(), p2_pixels.y(), rectangle_style);
    draw_layer_region_on(layer, plot_rectangle_cs.draw_object());
    return;
  }

  cs::Point<cs> const topleft_cs(plot_rectangle_cs.offset_x(), plot_rectangle_cs.offset_y());
  cs::Point<cs> const topright_cs(plot_rectangle_cs.offset_x() + plot_rectangle_cs.width(), plot_rectangle_cs.offset_y());
  cs::Point<cs> const bottomright_cs(plot_rectangle_cs.offset_x() + plot_rectangle_cs.width(), plot_rectangle_cs.offset_y() + plot_rectangle_cs.height());
  cs::Point<cs> const bottomleft_cs(plot_rectangle_cs.offset_x(), plot_rectangle_cs.offset_y() + plot_rectangle_cs.height());

  std::vector<cs::Point<csid::pixels>> points_pixels = {
    topleft_cs * cs_transform_pixels_,
    topright_cs * cs_transform_pixels_,
    bottomright_cs * cs_transform_pixels_,
    bottomleft_cs * cs_transform_pixels_
  };

  draw::PolylineStyle polyline_style(rectangle_style);
  plot_rectangle_cs.create_polyline_draw_object({}, std::move(points_pixels), polyline_style({.closed = true}));
  draw_layer_region_on(layer, plot_rectangle_cs.draw_object());
}

//----------------------------------------------------------------------------
// CircleHandle

template<CS cs>
void CoordinateMapper<cs>::add_circle(LayerPtr const& layer, draw::CircleStyle const& circle_style, CircleHandle const& plot_circle_cs)
{
  cs::Point<cs> const& center = plot_circle_cs.center();
  double const radius = plot_circle_cs.radius();

  cs::Point<csid::pixels> const center_pixels = center * cs_transform_pixels_;

  auto const [rx_x, rx_y] = cs_transform_pixels_.map_vector(radius, 0.0);
  auto const [ry_x, ry_y] = cs_transform_pixels_.map_vector(0.0, radius);
  double const radius_x_pixels = std::hypot(rx_x, rx_y);
  double const radius_y_pixels = std::hypot(ry_x, ry_y);

  cairowindow::Geometry geometry;
  switch (circle_style.position())
  {
    case draw::at_center:
      geometry = {center_pixels.x() - radius_x_pixels, center_pixels.y() - radius_y_pixels, 2.0 * radius_x_pixels, 2.0 * radius_y_pixels};
      break;
    case draw::at_corner:
    default:
      geometry = {center_pixels.x(), center_pixels.y(), radius_x_pixels, radius_y_pixels};
      break;
  }

  plot_circle_cs.create_draw_object({}, geometry, circle_style);
  draw_layer_region_on(layer, plot_circle_cs.draw_object());
}

//----------------------------------------------------------------------------
// ArcHandle

template<CS cs>
void CoordinateMapper<cs>::add_arc(LayerPtr const& layer, draw::ArcStyle const& arc_style, ArcHandle const& plot_arc_cs)
{
  cs::Point<cs> const& center_cs = plot_arc_cs.center();
  double const radius_cs = plot_arc_cs.radius();
  double const start_angle = plot_arc_cs.start_angle();
  double const end_angle = plot_arc_cs.end_angle();

  cs::Point<csid::pixels> const center_pixels = center_cs * cs_transform_pixels_;

  auto const [rx_x, rx_y] = cs_transform_pixels_.map_vector(radius_cs, 0.0);
  auto const [ry_x, ry_y] = cs_transform_pixels_.map_vector(0.0, radius_cs);
  double const radius_pixels = std::max(std::hypot(rx_x, rx_y), std::hypot(ry_x, ry_y));

  auto const [sx_x, sx_y] = cs_transform_pixels_.map_vector(std::cos(start_angle), std::sin(start_angle));
  auto const [ex_x, ex_y] = cs_transform_pixels_.map_vector(std::cos(end_angle), std::sin(end_angle));
  double const start_angle_pixels = std::atan2(sx_y, sx_x);
  double const end_angle_pixels = std::atan2(ex_y, ex_x);

  auto const [bx_x, bx_y] = cs_transform_pixels_.map_vector(1.0, 0.0);
  auto const [by_x, by_y] = cs_transform_pixels_.map_vector(0.0, 1.0);
  double const det = bx_x * by_y - bx_y * by_x;

  double pixel_start_angle = start_angle_pixels;
  double pixel_end_angle = end_angle_pixels;
  if (det < 0.0)
    std::swap(pixel_start_angle, pixel_end_angle);

  plot_arc_cs.create_draw_object({}, center_pixels.x(), center_pixels.y(), pixel_start_angle, pixel_end_angle, radius_pixels, arc_style);
  draw_layer_region_on(layer, plot_arc_cs.draw_object());
}

//----------------------------------------------------------------------------
// LineHandle (infinite)

template<CS cs>
void CoordinateMapper<cs>::add_clipped_line(LayerPtr const& layer, draw::LineStyle const& line_style, LineHandle const& plot_line_cs,
    math::Hyperblock<2> const& clip_rectangle_cs)
{
  cs::Direction<cs> const& direction = plot_line_cs.direction();
  cs::Point<cs> const& point = plot_line_cs.point();

  double normal_x = -direction.y();
  double normal_y = direction.x();
  math::Hyperplane<2> line({normal_x, normal_y}, -(normal_x * point.x() + normal_y * point.y()));
  auto intersections = clip_rectangle_cs.intersection_points(line);

  // Is the line outside the plot area?
  if (intersections.empty())
    return;

  constexpr math::Hyperblock<2>::IntersectionPointIndex first{size_t{0}};
  constexpr math::Hyperblock<2>::IntersectionPointIndex second{size_t{1}};

  double x1 = intersections[first].coordinate(0);
  double y1 = intersections[first].coordinate(1);
  double x2 = intersections[second].coordinate(0);
  double y2 = intersections[second].coordinate(1);

  cs::Point<cs> const intersection1_cs{x1, y1};
  cs::Point<cs> const intersection2_cs{x2, y2};

  cs::Point<csid::pixels> intersection1_pixels = intersection1_cs * cs_transform_pixels_;
  cs::Point<csid::pixels> intersection2_pixels = intersection2_cs * cs_transform_pixels_;

  plot_line_cs.create_draw_object({}, intersection1_pixels.x(), intersection1_pixels.y(), intersection2_pixels.x(), intersection2_pixels.y(), line_style);
  draw_layer_region_on(layer, plot_line_cs.draw_object());
}

//----------------------------------------------------------------------------
// LinePieceHandle

namespace detail {

inline void apply_line_extend(double& x1, double& y1, double& x2, double& y2, LineExtend line_extend, math::Hyperblock<2> const& clip_rectangle)
{
  if (line_extend == LineExtend::none)
    return;

  double dx = x2 - x1;
  double dy = y2 - y1;
  double normal_x = dy;
  double normal_y = -dx;
  math::Hyperplane<2> line({normal_x, normal_y}, -(normal_x * x1 + normal_y * y1));
  auto intersections = clip_rectangle.intersection_points(line);
  constexpr math::Hyperblock<2>::IntersectionPointIndex first{size_t{0}};
  constexpr math::Hyperblock<2>::IntersectionPointIndex second{size_t{1}};
  if (!intersections.empty())
  {
    // It is not known which intersection with the bounding rectangle ends up where in the intersections array.
    // Therefore look at the sign of the dot product between the line piece and the line between the two intersections.
    auto index_to = (dx * (intersections[second].coordinate(0) - intersections[first].coordinate(0)) +
                     dy * (intersections[second].coordinate(1) - intersections[first].coordinate(1))) < 0.0 ? second : first;
    if (line_extend == LineExtend::from || line_extend == LineExtend::both)
    {
      x1 = intersections[index_to].coordinate(0);
      y1 = intersections[index_to].coordinate(1);
    }
    if (line_extend == LineExtend::to || line_extend == LineExtend::both)
    {
      x2 = intersections[size_t{1} - index_to].coordinate(0);
      y2 = intersections[size_t{1} - index_to].coordinate(1);
    }
  }
}

} // namespace detail

template<CS cs>
void CoordinateMapper<cs>::add_clipped_line_piece(LayerPtr const& layer, draw::LineStyle const& line_style, LineExtend line_extend,
    LinePieceHandle const& plot_line_piece_cs, math::Hyperblock<2> const& clip_rectangle_cs)
{
  cs::Point<cs> const& from_cs = plot_line_piece_cs.from();
  cs::Point<cs> const& to_cs = plot_line_piece_cs.to();

  // apply_line_extend can change these initial values.
  double x1_cs = from_cs.x();
  double y1_cs = from_cs.y();
  double x2_cs = to_cs.x();
  double y2_cs = to_cs.y();
  detail::apply_line_extend(x1_cs, y1_cs, x2_cs, y2_cs, line_extend, clip_rectangle_cs);
  // Convert the result to pixels.
  cs::Point<csid::pixels> const p1_pixels = cs::Point<cs>{x1_cs, y1_cs} * cs_transform_pixels_;
  cs::Point<csid::pixels> const p2_pixels = cs::Point<cs>{x2_cs, y2_cs} * cs_transform_pixels_;

  plot_line_piece_cs.create_draw_object({}, p1_pixels.x(), p1_pixels.y(), p2_pixels.x(), p2_pixels.y(), line_style);
  draw_layer_region_on(layer, plot_line_piece_cs.draw_object());
}

//----------------------------------------------------------------------------
// ConnectorHandle

template<CS cs>
void CoordinateMapper<cs>::add_connector(LayerPtr const& layer, draw::ConnectorStyle const& connector_style, ConnectorHandle const& plot_connector_cs)
{
  cs::Point<cs> const& from = plot_connector_cs.from();
  cs::Point<cs> const& to = plot_connector_cs.to();

  auto const arrow_head_shape_from = static_cast<draw::Connector::ArrowHeadShape>(plot_connector_cs.arrow_head_shape_from());
  auto const arrow_head_shape_to = static_cast<draw::Connector::ArrowHeadShape>(plot_connector_cs.arrow_head_shape_to());

  cs::Point<csid::pixels> const from_pixels = from * cs_transform_pixels_;
  cs::Point<csid::pixels> const to_pixels = to * cs_transform_pixels_;

  plot_connector_cs.create_draw_object({},
      from_pixels.x(), from_pixels.y(), to_pixels.x(), to_pixels.y(),
      connector_style, arrow_head_shape_from, arrow_head_shape_to);

  if (need_print_)
    layer->start_printing_to(svg_cr_);
  layer->draw(plot_connector_cs.draw_object());
  plot_connector_cs.draw_object()->draw_arrow_heads(layer);
  if (need_print_)
    layer->stop_printing();
}

} // namespace cairowindow
