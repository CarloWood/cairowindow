#pragma once

#include "Layer.h"
#include "plot/Rectangle.h"
#include "plot/Line.h"
#include "cs/Vector.h"
#include "cs/Range.h"
#include "cs/NiceDelta.h"
#include "draw/Point.h"
#include "draw/Line.h"
#include "draw/Rectangle.h"
#include "draw/Polyline.h"
#include "draw/Text.h"
#include "draw/Arc.h"
#include "draw/Circle.h"
#include "draw/Connector.h"
#include "intersection_points.h"
#include "CoordinateMapper.h"
#include "math/Hyperblock.h"
#include "math/Transform.h"
#include <boost/intrusive_ptr.hpp>
#include <cmath>
#include <string>
#include <vector>
#include "debug.h"

namespace cairowindow {

//using PointStyle = cwin::draw::PointStyle;
//using ConnectorStyle = cwin::draw::ConnectorStyle;
//using RectangleStyle = cwin::draw::RectangleStyle;
//using CircleStyle = cwin::draw::CircleStyle;
//using ArcStyle = cwin::draw::ArcStyle;
//using TextStyle = cwin::draw::TextStyle;

//---------------------------------------------------------------------------
// CoordinateSystem
//
// Draw the 2D coordinate system `cs` with axes, tick marks and labels from a
// given math::Transform<cs, csid::pixels>.
//
// Usage
// -----
//
//   CoordinateSystem<csid::foo> foo_cs(cs_transform_pixels, window.geometry());
//   foo_cs.display(layer, axis_style);
//
// where `cs_transform_pixels` is a math::Transform<csid::foo, csid::pixels> that converts
// from a logical coordinate system `foo` to window pixels.
//
// Lifetime and ownership
// ----------------------
//
// The CoordinateSystem object must be kept alive for as long as you
// want the axes and tick labels to remain visible.
//
// API overview
// ------------
// - Constructor
//     CoordinateSystem(math::Transform<cs, csid::pixels> cs_transform_pixels, Geometry window_geometry);
//
//   `cs_transform_pixels` converts from your logical coordinates `cs`
//   to window pixels. The constructor uses this transform to compute
//   where the x and y axes intersect the window rectangle
//   `window_geometry`. The visible part of each
//   axis is stored and `set_range` is called with the corresponding
//   range in `cs` coordinates. Tick spacing and label formatting are
//   chosen automatically using NiceDelta<cs>.
//
// - void display(LayerPtr const& layer, LineStyle const& axis_style);
//     Draws the axes, tick marks and labels on `layer`.
//     Call once per CoordinateSystem instance.
//
template<CS cs>
class CoordinateSystem : public CoordinateMapper<cs>
{
  static constexpr int x_axis = 0;
  static constexpr int y_axis = 1;
  static constexpr int number_of_axes = 2;

 private:
  cs::Point<csid::pixels> csOrigin_pixels_;                                     // The origin in pixels.
  std::array<cs::Direction<csid::pixels>, number_of_axes> csAxisDirection_;     // The direction of the x-axis and y-axis (in csid::pixels).
  std::vector<std::shared_ptr<draw::Line>> lines_;                              // To keep drawn lines alive.
  std::vector<std::shared_ptr<draw::Text>> texts_;                              // To keep drawn texts alive.
  std::array<cs::LinePiece<csid::pixels>, number_of_axes> line_piece_;          // The visible part of the axes (in csid::pixels).
  Geometry window_geometry_;                                                    // Geometry used to calculate visible axis segments.

 private:
  using LayerPtr = boost::intrusive_ptr<Layer>;

  using PointHandle     = plot::cs::Point<cs>;
  using RectangleHandle = plot::cs::Rectangle<cs>;
  using LineHandle      = plot::cs::Line<cs>;
  using LinePieceHandle = plot::cs::LinePiece<cs>;

 protected:
  // Bring base class members into scope (dependent base).
  using CoordinateMapper<cs>::cs_transform_pixels_;

 public:
//  using LinePiecePtr = std::shared_ptr<LinePiece<cs>>;
//  using CirclePtr = std::shared_ptr<Circle<cs>>;
//  using ArcPtr = std::shared_ptr<Arc<cs>>;
//  using TextPtr = std::shared_ptr<Text<cs>>;

/*  std::shared_ptr<Text> xlabel_;
  std::shared_ptr<Text> ylabel_;*/
  std::array<cs::Range<cs>, number_of_axes> range_{
    {{0.0, 0.0}, {0.0, 0.0}}                                                    // Zero means: not visible.
  };
  std::array<cs::NiceDelta<cs>, number_of_axes> range_ticks_;                   // The number of tick marks on the visible segment of the respective axis.
                                                                                // Invalid (default constructed) means: don't draw ticks.
/*  std::array<std::vector<std::shared_ptr<Text>>, number_of_axes> labels_;*/

 public:
  CoordinateSystem(math::Transform<cs, csid::pixels> const reference_transform, Geometry window_geometry);

  ~CoordinateSystem()
  {
    DoutEntering(dc::notice, "CoordinateSystem::~CoordinateSystem() [" << this << "]");
  }

  Geometry const& geometry() const override
  {
    return window_geometry_;
  }

  void set_range(int axis, cs::Range<cs> range)
  {
    DoutEntering(dc::notice, "CoordinateSystem::set_range(" << axis << ", " << range << ") [" << this << "]");
    range_[axis] = range;
    range_ticks_[axis] = cs::NiceDelta<cs>{range};
    Dout(dc::notice, "range_[" << axis << "] = " << range_[axis] << "; range_ticks_[" << axis << "] = " << range_ticks_[axis]);
  }

  // Accessors.
  Range const& xrange() const { return range_[x_axis]; }
  Range const& yrange() const { return range_[y_axis]; }

#if 0
  double convert_x(double x) const;
  double convert_y(double y) const;
  Pixel convert_to_pixel(Point const& point) const;

  double convert_from_pixel_x(double pixel_x) const;
  double convert_from_pixel_y(double pixel_y) const;
  Point convert_from_pixel(Pixel const& pixel) const;

  double convert_horizontal_offset_from_pixel(double pixel_offset_x) const;
  double convert_vertical_offset_from_pixel(double pixel_offset_y) const;

  void convert_to_pixels(Point const* data_in, Pixel* data_out, std::size_t size);

  template<std::size_t size>
  [[gnu::always_inline]] void convert_to_pixels(std::array<Point, size> const& in, std::array<Pixel, size>& out)
  {
    convert_to_pixels(in.data(), out.data(), size);
  }
#endif

  //--------------------------------------------------------------------------
  // Point

 public:
  // Create and draw a point on layer at x,y using point_style.
  [[nodiscard]] PointHandle create_point(LayerPtr const& layer, draw::PointStyle const& point_style, cs::Point<cs> const& point_cs)
  {
    PointHandle plot_point_cs{point_cs};
    add_point(layer, point_style, plot_point_cs);
    return plot_point_cs;
  }

#if 0
  //--------------------------------------------------------------------------
  // Connector

  // Add and draw cs_connector on layer, using line_style, fill_color, arrow_head_shape_from and arrow_head_shape_to.
  void add_connector(LayerPtr const& layer, ConnectorStyle const& style, ConnectorPtr const& cs_connector);

 private:
  void add_connector(LayerPtr const& layer, ConnectorStyle const& style, ConnectorPtr&& cs_connector);

 private:
  template<typename... Args>
  [[nodiscard]] ConnectorPtr create_connector_helper(
      Connector::ArrowHeadShape& default_arrow_head_shape_from,
      Connector::ArrowHeadShape& default_arrow_head_shape_to,
      LayerPtr const& layer,
      ConnectorStyle const& connector_style, std::tuple<Args...>&& args)
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
      ConnectorPtr cs_connector = std::apply([&](auto&&... unpacked_args) -> ConnectorPtr {
        return std::make_shared<Connector<cs>>(std::forward<decltype(unpacked_args)>(unpacked_args)..., default_arrow_head_shape_from, default_arrow_head_shape_to);
      }, std::move(args));

      add_connector(layer, connector_style, cs_connector);
      return cs_connector;
    }
  }

 public:
  // Create and draw a connector on layer, using args... and line_style.
  // Args can optionally start with zero, one or two ArrowHeadShape arguments.
  template<typename... Args>
  [[nodiscard]] ConnectorPtr create_connector(LayerPtr const& layer,
      ConnectorStyle const& connector_style, Args&&... args)
  {
    Connector::ArrowHeadShape arrow_head_shape_from = Connector::no_arrow;
    Connector::ArrowHeadShape arrow_head_shape_to = Connector::open_arrow;
    return create_connector_helper(arrow_head_shape_from, arrow_head_shape_to, layer,
        connector_style, std::make_tuple(std::forward<Args>(args)...));
  }
#endif

  //--------------------------------------------------------------------------
  // Line

  // Add and draw plot_line using line_style.
  void add_line(LayerPtr const& layer, draw::LineStyle const& line_style, LineHandle const& plot_line_cs);

 public:
  // Create and draw a line through point in direction using line_style.
  template<typename... Args>
  [[nodiscard]] LineHandle create_line(LayerPtr const& layer, draw::LineStyle const& line_style, Args&&... args)
  {
    LineHandle plot_line_cs(std::forward<Args>(args)...);
    add_line(layer, line_style, plot_line_cs);
    return plot_line_cs;
  }

 public:
  // Create and draw a rectangle on layer, using args... and rectangle_style.
  template<typename... Args>
  [[nodiscard]] RectangleHandle create_rectangle(LayerPtr const& layer, draw::RectangleStyle const& rectangle_style, Args&&... args)
  {
    RectangleHandle plot_rectangle_cs(std::forward<Args>(args)...);
    add_rectangle(layer, rectangle_style, plot_rectangle_cs);
    return plot_rectangle_cs;
  }

  //--------------------------------------------------------------------------
  // LinePiece

  // Add and draw plot_line_piece_cs on layer using line_style and line_extend.
  void add_line_piece(LayerPtr const& layer, draw::LineStyle const& line_style, LineExtend line_extend, LinePieceHandle const& plot_line_piece_cs);

 public:
  // Create and draw a line piece between points from and to using line_style and line_extend.
  template<typename... Args>
  [[nodiscard]] LinePieceHandle create_line_piece(LayerPtr const& layer, draw::LineStyle const& line_style, LineExtend line_extend, Args&&... args)
  {
    LinePieceHandle plot_line_piece_cs(std::forward<Args>(args)...);
    add_line_piece(layer, line_style, line_extend, plot_line_piece_cs);
    return plot_line_piece_cs;
  }

  // Same, but without a line_extend.
  template<typename... Args>
  [[nodiscard]] LinePieceHandle create_line_piece(LayerPtr const& layer, draw::LineStyle const& line_style, Args&&... args)
  {
    LinePieceHandle plot_line_piece_cs(std::forward<Args>(args)...);
    add_line_piece(layer, line_style, LineExtend::none, plot_line_piece_cs);
    return plot_line_piece_cs;
  }

#if 0
  //--------------------------------------------------------------------------
  // Circle

  // Add and draw cs_circle on layer with center and radius using circle_style.
  void add_circle(LayerPtr const& layer, CircleStyle const& circle_style, CirclePtr const& cs_circle);

 private:
  void add_circle(LayerPtr const& layer, CircleStyle const& circle_style, CirclePtr&& cs_circle);

 public:
  // Create and draw a circle on layer with center and radius using circle_style.
  template<typename... Args>
  [[nodiscard]] CirclePtr create_circle(LayerPtr const& layer, CircleStyle const& circle_style, Args&&... args)
  {
    CirclePtr cs_circle = std::make_shared<Circle<cs>>(std::forward<Args>(args)...);
    add_circle(layer, circle_style, cs_circle);
    return cs_circle;
  }

  // Same as above but use line_style (no fill_color).
  template<typename... Args>
  [[nodiscard]] CirclePtr create_circle(LayerPtr const& layer, LineStyle const& line_style, Args&&... args)
  {
    return create_circle(layer, CircleStyle({.line_color = line_style.line_color(), .line_width = line_style.line_width()}), std::forward<Args>(args)...);
  }
#endif

#if 0
  //--------------------------------------------------------------------------
  // Arc

  // Add and draw cs_arc on layer using arc_style.
  void add_arc(LayerPtr const& layer, ArcStyle const& arc_style, ArcPtr const& cs_arc);

 private:
  void add_arc(LayerPtr const& layer, ArcStyle const& arc_style, ArcPtr&& cs_arc);

 public:
  // Create and draw an arc on layer width center, radius and start- and end_angle, using arc_style.
  template<typename... Args>
  [[nodiscard]] ArcPtr create_arc(LayerPtr const& layer, ArcStyle const& arc_style, Args&&... args)
  {
    ArcPtr cs_arc = std::make_shared<Arc<cs>>(std::forward<Args>(args)...);
    add_arc(layer, arc_style, cs_arc);
    return cs_arc;
  }
#endif

#if 0
  //--------------------------------------------------------------------------
  // Text

  // Add and draw cs_text on layer using text_style.
  void add_text(LayerPtr const& layer, TextStyle const& text_style, TextPtr const& cs_text);

 private:
  void add_text(LayerPtr const& layer, TextStyle const& text_style, TextPtr&& cs_text);

 public:
  // Create and draw text on layer at position using text_style.
  [[nodiscard]] TextPtr create_text(LayerPtr const& layer, TextStyle const& text_style, cs::Point<cs> const& position, std::string const& text)
  {
    TextPtr cs_text = std::make_shared<Text<cs>>(convert_to_pixel(position), text);
    add_text(layer, text_style, cs_text);
    return cs_text;
  }

  // Same, but using pixel coordinates.
  [[nodiscard]] TextPtr create_text(LayerPtr const& layer, TextStyle const& text_style, cs::Point<pixel> position, std::string const& text)
  {
    TextPtr cs_text = std::make_shared<Text<cs>>(position, text);
    add_text(layer, text_style, cs_text);
    return cs_text;
  }
#endif

  void display(LayerPtr const& layer, draw::LineStyle const& axis_style);

 private:
//  void apply_line_extend(double& x1, double& y1, double& x2, double& y2, LineExtend line_extend);
};

namespace detail {

// Calculate the intersection between a line and a rectangle and return
// the number of intersection points found (0 or 2) and an array containing
// those intersection points, if any.
//
// The line has a direction; this direction will correspond with the
// direction from the point returned as index 0 of the array, to the
// point returned as index 1.
//
// The points are returned as a cs::Point<cs>. The caller is responsible to
// make sure that the rectangle uses that same coordinate system.
template<CS cs>
std::tuple<int, std::array<cs::Point<cs>, 2>> intersect(cs::Line<cs> line_cs, math::Hyperblock<2> rectangle_cs)
{
//  DoutEntering(dc::notice, "detail::intersect(" << line_cs << ", " << rectangle_cs << ")");

  double normal_x = -line_cs.direction().y();
  double normal_y = line_cs.direction().x();
  math::Hyperplane<2> line({normal_x, normal_y}, -(normal_x * line_cs.point().x() + normal_y * line_cs.point().y()));
  auto intersections_cs = rectangle_cs.intersection_points(line);

  // Is the line outside the window?
  if (intersections_cs.empty())
    return {0, {}};

  ASSERT(intersections_cs.size() == 2);

  constexpr math::Hyperblock<2>::IntersectionPointIndex first{size_t{0}};
  constexpr math::Hyperblock<2>::IntersectionPointIndex second{size_t{1}};

  // Return the two points where the line_cs intersects with the rectangle_cs.
  return {
    2, {
      cs::Point<cs>(intersections_cs[first].coordinate(0), intersections_cs[first].coordinate(1)),
      cs::Point<cs>(intersections_cs[second].coordinate(0), intersections_cs[second].coordinate(1))
    }
  };
}

} // namespace detail

template<CS cs>
CoordinateSystem<cs>::CoordinateSystem(math::Transform<cs, csid::pixels> const cs_transform_pixels, Geometry window_geometry) :
  CoordinateMapper<cs>(cs_transform_pixels), window_geometry_(window_geometry)
{
  DoutEntering(dc::notice, "CoordinateSystem::CoordinateSystem(" << cs_transform_pixels << ") [" << this << "]");

  // Calculate where the cs-axis intersect with the window geometry.

  // Construct three points on the CS axis.
  cs::Point<cs> csOrigin_cs;
  cs::Point<cs> csXAxisUnit_cs(1, 0);
  cs::Point<cs> csYAxisUnit_cs(0, 1);
  // Convert them to pixel coordinates.
  csOrigin_pixels_ = csOrigin_cs * cs_transform_pixels;
  cs::Point<csid::pixels> csXAxisUnit_pixels = csXAxisUnit_cs * cs_transform_pixels;
  cs::Point<csid::pixels> csYAxisUnit_pixels = csYAxisUnit_cs * cs_transform_pixels;
  // Store the directions (in pixels coordinates).
  csAxisDirection_[x_axis] = cs::Direction<csid::pixels>(csOrigin_pixels_, csXAxisUnit_pixels);
  csAxisDirection_[y_axis] = cs::Direction<csid::pixels>(csOrigin_pixels_, csYAxisUnit_pixels);

  // The inverse transform.
  auto const& pixels_transform_cs = cs_transform_pixels.inverse();

  // Calcuate the visible length of each CS axis, in pixels.
  for (int axis = x_axis; axis <= y_axis; ++axis)
  {
    // Determine where the axis intersects with the window rectangle (everything in pixels).
    auto [number_of_intersection_points, intersection_point_pixels] = detail::intersect<csid::pixels>(
        {csOrigin_pixels_, csAxisDirection_[axis]},     // The axis (pointing in the direction csAxisDirection_).
        {{window_geometry_.offset_x(), window_geometry_.offset_y()}, {window_geometry_.width(), window_geometry_.height()}});   // The window rectangle.

    // Is the line outside the window?
    if (number_of_intersection_points < 2)
    {
      cs::Point<csid::pixels> const origin(0, 0);
      line_piece_[axis] = cs::LinePiece<csid::pixels>{origin, origin};  // Use twice the same point to encode that this axis is not visible within the window.
      continue;
    }

    // See detail::intersect.
    constexpr int from = 0;     // The index into intersection_point_pixels where the negative side of the axis intersects with the window rectangle.
    constexpr int to = 1;       // Same, but the positive side of the axis.

    line_piece_[axis] = cs::LinePiece<csid::pixels>(intersection_point_pixels[from], intersection_point_pixels[to]);
    Dout(dc::notice, "line_piece_[" << axis << "] = " << line_piece_[axis]);

    // Convert the intersection points back to cs.
    cs::Point<cs> const from_cs = intersection_point_pixels[from] * pixels_transform_cs;
    cs::Point<cs> const   to_cs =   intersection_point_pixels[to] * pixels_transform_cs;
    Dout(dc::notice, "from_cs = " << from_cs);
    Dout(dc::notice, "to_cs = " << to_cs);
    // Extract the minimum and maximum values of the visible range.
    double min = (axis == x_axis) ? from_cs.x() : from_cs.y();
    double max = (axis == x_axis) ?   to_cs.x() :   to_cs.y();
    ASSERT(min < max);
    // Set the range on each (visible) axis.
    set_range(axis, {min, max});
  }
}

template<CS cs>
void CoordinateSystem<cs>::display(LayerPtr const& layer, draw::LineStyle const& axis_style)
{
  DoutEntering(dc::notice, "CoordinateSystem<" << utils::to_string(cs) << ">::display(layer, axis_style)");

  using namespace cairowindow;

  ASSERT(texts_.empty());
  for (int axis = x_axis; axis <= y_axis; ++axis)
  {
    // If the range is empty then min = max = 0 and size() will return zero exactly.
    if (range_[axis].size() == 0.0)     // Not visible?
      continue;
    // Draw the piece of the axis that is visible.
    lines_.emplace_back(std::make_shared<draw::Line>(
          line_piece_[axis].from().x(), line_piece_[axis].from().y(),
          line_piece_[axis].to().x(), line_piece_[axis].to().y(),
          axis_style));
    layer->draw(lines_.back());

    if (range_ticks_[axis].is_invalid())
      continue;

    // Draw the tick marks.
    double const delta_cs = range_ticks_[axis].value();
    int k_min = std::ceil(range_[axis].min() / delta_cs);
    int k_max = std::floor(range_[axis].max() / delta_cs);
    auto const axis_direction = csAxisDirection_[axis];
    bool const axis_prefers_parallel = std::abs(axis_direction.x()) >= std::abs(axis_direction.y());
    auto const axis_angle = axis_direction.as_angle();
    double const pi = std::acos(-1.0);
    auto normalize_readable = [pi](double angle) {
      if (angle > 0.5 * pi)
        angle -= pi;
      else if (angle <= -0.5 * pi)
        angle += pi;
      return angle;
    };
    for (int k = k_min; k <= k_max; ++k)
    {
      if (k == 0)
        continue;
      double value_cs = k * delta_cs;
      cs::Point<cs> tick_cs{axis == x_axis ? value_cs : 0.0, axis == y_axis ? value_cs : 0.0};
      cs::Point<csid::pixels> tick_pixels = tick_cs * cs_transform_pixels_;
      // Unit vector pointing into the positive direction of axis.
      cs::Direction<csid::pixels> axis_pixels{csOrigin_pixels_, tick_pixels};
      cs::Direction<csid::pixels> axis_tickmark_pixels =
          (axis == x_axis) == (k < 0) ? axis_pixels.rotated_90_degrees() : axis_pixels.rotated_270_degrees();
      cs::Point<csid::pixels> tick_end_pixels = tick_pixels + cs::Vector<csid::pixels>{axis_tickmark_pixels, 5.0};
      lines_.emplace_back(std::make_shared<draw::Line>(
            tick_pixels.x(), tick_pixels.y(),
            tick_end_pixels.x(), tick_end_pixels.y(),
            axis_style));
      layer->draw(lines_.back());

      cs::Point<csid::pixels> text_anchor_pixels = tick_pixels + cs::Vector<csid::pixels>{axis_tickmark_pixels, 10.0};

      std::string label = range_ticks_[axis].label(k);

      double rotation = axis_angle;
      draw::TextPosition position;
      if (axis_prefers_parallel)
        position = axis_tickmark_pixels.y() < 0 ? draw::centered_above : draw::centered_below;
      else
      {
        rotation += 0.5 * pi;
        position = axis_tickmark_pixels.x() < 0 ? draw::centered_left_of : draw::centered_right_of;
      }
      rotation = normalize_readable(rotation);

      draw::TextStyle text_style({
          .position = position,
          .color = axis_style.line_color(),
          .rotation = rotation
      });

      texts_.emplace_back(std::make_shared<draw::Text>(
            label,
            text_anchor_pixels.x(),
            text_anchor_pixels.y(),
            text_style));
      layer->draw(texts_.back());
    }
  }
}

//--------------------------------------------------------------------------
// Line

template<CS cs>
void CoordinateSystem<cs>::add_line(LayerPtr const& layer, draw::LineStyle const& line_style, LineHandle const& plot_line_cs)
{
  cs::Direction<cs> const& direction_cs = plot_line_cs.direction();
  cs::Point<cs> const& point_cs = plot_line_cs.point();

  // Convert the line through `point` with direction `direction` into a HyperPlane.
  double normal_x = -direction_cs.y();
  double normal_y = direction_cs.x();
  math::Hyperplane<2> line_cs({normal_x, normal_y}, -(normal_x * point_cs.x()+ normal_y * point_cs.y()));

  // Construct a Rectangle aligned with the axes that captures the range of the axes.
  // Make the rectangle twice as big relative to its center.
  cs::Rectangle<cs> rectangle_cs{
    range_[x_axis].center() - range_[x_axis].size(),
    range_[y_axis].center() - range_[y_axis].size(),
    2 * range_[x_axis].size(),
    2 * range_[y_axis].size()
  };
  math::Hyperblock<2> const clip_rectangle_cs(
      {rectangle_cs.offset_x(), rectangle_cs.offset_y()},
      {rectangle_cs.offset_x() + rectangle_cs.width(), rectangle_cs.offset_y() + rectangle_cs.height()});

  this->add_clipped_line(layer, line_style, plot_line_cs, clip_rectangle_cs);
}

//--------------------------------------------------------------------------
// LinePiece

template<CS cs>
void CoordinateSystem<cs>::add_line_piece(LayerPtr const& layer, draw::LineStyle const& line_style, LineExtend line_extend, LinePieceHandle const& plot_line_piece_cs)
{
  cs::Point<cs> const& from_cs = plot_line_piece_cs.from();
  cs::Point<cs> const& to_cs = plot_line_piece_cs.to();
  cs::Direction<cs> const direction_cs(from_cs, to_cs);

  // Convert the line through `from_cs` with direction `direction` into a HyperPlane.
  double normal_x = -direction_cs.y();
  double normal_y = direction_cs.x();
  math::Hyperplane<2> line_cs({normal_x, normal_y}, -(normal_x * from_cs.x()+ normal_y * from_cs.y()));

  // Construct a Rectangle aligned with the axes that captures the range of the axes.
  // Make the rectangle twice as big relative to its center.
  cs::Rectangle<cs> rectangle_cs{
    range_[x_axis].center() - range_[x_axis].size(),
    range_[y_axis].center() - range_[y_axis].size(),
    2 * range_[x_axis].size(),
    2 * range_[y_axis].size()
  };
  math::Hyperblock<2> const clip_rectangle_cs(
      {rectangle_cs.offset_x(), rectangle_cs.offset_y()},
      {rectangle_cs.offset_x() + rectangle_cs.width(), rectangle_cs.offset_y() + rectangle_cs.height()});

  this->add_clipped_line_piece(layer, line_style, line_extend, plot_line_piece_cs, clip_rectangle_cs);
}

#if 0
//--------------------------------------------------------------------------
// Text

template<CS cs>
void CoordinateSystem<cs>::add_text(LayerPtr const& layer, TextStyle const& text_style, TextPtr const& cs_text)
{
  Pixel position = plot_text.position();
  std::string const& text = plot_text.text();

  plot_text.draw_object_ = std::make_shared<draw::Text>(text, position.x(), position.y(), text_style);
  draw_layer_region_on(layer, plot_text.draw_object_);
}
#endif

} // namespace cairowindow
