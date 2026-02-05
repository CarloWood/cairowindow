#pragma once

#include "Printable.h"
#include "cs/Point.h"
#include "cs/TransformOperators.h"
#include "plot/Line.h"
#include "plot/Point.h"
#include "plot/Rectangle.h"
#include "draw/Rectangle.h"
#include "draw/Polyline.h"
#include "math/Hyperblock.h"
#include "math/Transform.h"
#include <functional>

namespace cairowindow {

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

  // Add and draw cs_point on layer using point_style.
  void add_point(LayerPtr const& layer, draw::PointStyle const& point_style, PointHandle const& plot_point_cs);

  // Called by Window::update_grabbed through the lambda defined in Window::register_draggable<cs> when a Draggable plot_point_cs was moved to (pixel_x, pixel_y).
  Geometry update_grabbed(plot::cs::Point<cs>* plot_point_cs, double pixel_x, double pixel_y,
      std::function<cs::Point<cs> (cs::Point<cs> const&)> const& restriction);

  //--------------------------------------------------------------------------
  // Rectangle

  // Add and draw cs_rectangle on layer using rectangle_style. Handles rotation by drawing a polyline.
  void add_rectangle(LayerPtr const& layer, draw::RectangleStyle const& rectangle_style, RectangleHandle const& plot_rectangle_cs);

  //--------------------------------------------------------------------------
  // Line (infinite, clipped)

  // Add and draw cs_line on layer using line_style, clipped to clip_rectangle_cs.
  void add_clipped_line(LayerPtr const& layer, draw::LineStyle const& line_style, LineHandle const& plot_line_cs, math::Hyperblock<2> const& clip_rectangle_cs);
};

//----------------------------------------------------------------------------
// CoordinateMapper::Point

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
// CoordinateMapper::Rectangle

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
// CoordinateMapper::Line (infinite)

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

} // namespace cairowindow
