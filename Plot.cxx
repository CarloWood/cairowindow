#include "sys.h"
#include "Plot.h"
#include "Layer.h"
#include "math/Hyperblock.h"
#include "utils/square.h"
#include <limits>
#include <cmath>
#include <iomanip>

namespace cairowindow::plot {

// Calculate the axes geometry from the full plot geometry.
cairowindow::Geometry Plot::axes_geometry(cairowindow::Geometry const& geometry, double axes_line_width)
{
  // Use fixed-size margins for now.
  constexpr int top_margin = 50;
  constexpr int bottom_margin = 60;
  constexpr int left_margin = 100;
  constexpr int right_margin = 20;

  double top_left_x = geometry.offset_x() + left_margin;
  double top_left_y = geometry.offset_y() + top_margin;
  double bottom_right_x = geometry.offset_x() + geometry.width() - right_margin;
  double bottom_right_y = geometry.offset_y() + geometry.height() - bottom_margin;

  int lwi = std::round(axes_line_width);
  if (lwi == axes_line_width)
  {
    if ((lwi & 1) == 0)
    {
      top_left_x = std::round(top_left_x);
      top_left_y = std::round(top_left_y);
      bottom_right_x = std::round(bottom_right_x);
      bottom_right_y = std::round(bottom_right_y);
    }
    else
    {
      top_left_x = std::ceil(top_left_x) - 0.5;
      top_left_y = std::ceil(top_left_y) - 0.5;
      bottom_right_x = std::floor(bottom_right_x) + 0.5;
      bottom_right_y = std::floor(bottom_right_y) + 0.5;
    }
  }

  return { top_left_x, top_left_y, bottom_right_x - top_left_x, bottom_right_y - top_left_y };
}

void Plot::add_to(boost::intrusive_ptr<Layer> const& layer, bool keep_ratio)
{
  if (keep_ratio)
  {
    // Fix plot_area_ geometry.
    cairowindow::Geometry const& geometry = plot_area_.geometry();
    double pixels_per_x_unit = geometry.width() / range_[x_axis].size();
    double pixels_per_y_unit = geometry.height() / range_[y_axis].size();
    double required_scale = pixels_per_y_unit / pixels_per_x_unit;
    if (required_scale < 1.0)
    {
      if (xlabel_)
        xlabel_->rel_move_to(-0.5 * geometry.width() * (1.0 - required_scale), 0.0);
      plot_area_.set_geometry({ geometry.offset_x(), geometry.offset_y(), geometry.width() * required_scale, geometry.height() });
    }
    else
    {
      if (ylabel_)
        ylabel_->rel_move_to(0.0, -0.5 * geometry.height() * (1.0 - 1.0 / required_scale));
      plot_area_.set_geometry({ geometry.offset_x(), geometry.offset_y(), geometry.width(), geometry.height() / required_scale });
    }
  }

  // Set a title.
  if (title_)
  {
    title_->move_to(plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_->style().offset());
    draw_layer_region_on(layer, title_);
  }

  // Set ranges on the plot area and draw it.
  for (int axis = x_axis; axis <= y_axis; ++axis)
    plot_area_.set_range(axis, range_[axis].min(), range_[axis].max(), range_ticks_[axis]);
  draw_multi_region_on(layer, &plot_area_);

  // Draw axis labels.
  double max_height = 0;
  double max_width = 0;
  for (int axis = x_axis; axis <= y_axis; ++axis)
  {
    double label = range_[axis].min();
    double delta = range_[axis].size() / range_ticks_[axis];
    int precision = -std::floor(std::log10(delta) + 1e-6);
    double x = plot_area_.geometry().offset_x();
    double y = plot_area_.geometry().offset_y() + plot_area_.geometry().height();
    if (axis == x_axis)
      y += XLabelStyleDefaults::offset;
    else
      x -= XLabelStyleDefaults::offset;
    for (int tick = 0; tick <= range_ticks_[axis]; ++tick)
    {
      std::ostringstream label_str;
      if (precision > 0)
        label_str << std::fixed << std::setprecision(precision);
      label_str << label;
      labels_[axis].emplace_back(std::make_shared<draw::Text>(label_str.str(),
            x, y, draw::LabelStyle({.position = (axis == x_axis) ? draw::centered_below : draw::centered_left_of})));
      draw_layer_region_on(layer, labels_[axis].back());
      StrokeExtents label_extents = labels_[axis].back()->stroke_extents();
      if (axis == x_axis)
        max_height = std::max(max_height, label_extents.height());
      else
        max_width = std::max(max_width, label_extents.width());

      label += delta;
      if (axis == x_axis)
        x += plot_area_.geometry().width() / range_ticks_[axis];
      else
        y -= plot_area_.geometry().height() / range_ticks_[axis];
    }
  }

  if (xlabel_)
  {
    xlabel_->rel_move_to(0, max_height);
    draw_layer_region_on(layer, xlabel_);
  }

  if (ylabel_)
  {
    ylabel_->rel_move_to(-max_width, 0);
    draw_layer_region_on(layer, ylabel_);
  }

  // Register this plot and its geometry with the associated Window so that we can find which printable is under the mouse if needed.
  layer->window()->add_printable(this);
}

double Plot::convert_x(double x) const
{
  cairowindow::Geometry const& g = plot_area_.geometry();
  x = (x - range_[x_axis].min()) / range_[x_axis].size();
  x *= g.width();
  x += g.offset_x();
  return x;
}

double Plot::convert_y(double y) const
{
  cairowindow::Geometry const& g = plot_area_.geometry();
  y = (range_[y_axis].max() - y) / range_[y_axis].size();
  y *= g.height();
  y += g.offset_y();
  return y;
}

Pixel Plot::convert_to_pixel(cairowindow::Point const& point) const
{
  double x = point.x();
  double y = point.y();

  cairowindow::Geometry const& g = plot_area_.geometry();

  x -= range_[x_axis].min();
  x /= range_[x_axis].size();
  x *= g.width();
  x += g.offset_x();

  y = range_[y_axis].max() - y;
  y /= range_[y_axis].size();
  y *= g.height();
  y += g.offset_y();

  return Pixel{x, y};
}

double Plot::convert_from_pixel_x(double pixel_x) const
{
  double x = pixel_x;

  cairowindow::Geometry const& g = plot_area_.geometry();

  x -= g.offset_x();
  x /= g.width();
  x *= range_[x_axis].size();
  x += range_[x_axis].min();

  return x;
}

double Plot::convert_from_pixel_y(double pixel_y) const
{
  double y = pixel_y;

  cairowindow::Geometry const& g = plot_area_.geometry();

  y -= g.offset_y();
  y /= g.height();
  y *= range_[y_axis].size();
  y = range_[y_axis].max() - y;

  return y;
}

cairowindow::Point Plot::convert_from_pixel(Pixel const& pixel) const
{
  double x = pixel.x();
  double y = pixel.y();

  cairowindow::Geometry const& g = plot_area_.geometry();

  x -= g.offset_x();
  x /= g.width();
  x *= range_[x_axis].size();
  x += range_[x_axis].min();

  y -= g.offset_y();
  y /= g.height();
  y *= range_[y_axis].size();
  y = range_[y_axis].max() - y;

  return {x, y};
}

double Plot::convert_horizontal_offset_from_pixel(double pixel_offset_x) const
{
  cairowindow::Geometry const& g = plot_area_.geometry();
  return pixel_offset_x / g.width() * range_[x_axis].size();
}

double Plot::convert_vertical_offset_from_pixel(double pixel_offset_y) const
{
  cairowindow::Geometry const& g = plot_area_.geometry();
  return pixel_offset_y / g.height() * range_[y_axis].size();
}

void Plot::convert_to_pixels(cairowindow::Point const* data_in, Pixel* data_out, std::size_t size)
{
  cairowindow::Geometry const& g = plot_area_.geometry();
  double const x_offset = -range_[x_axis].min();
  double const y_offset = -range_[y_axis].max();
  double const x_scale = g.width() / range_[x_axis].size();
  double const y_scale = -g.height() / range_[y_axis].size();

  for (std::size_t i = 0; i < size; ++i)
  {
    double x = data_in[i].x();
    double y = data_in[i].y();

    x += x_offset;
    y += y_offset;

    x *= x_scale;
    y *= y_scale;

    x += g.offset_x();
    y += g.offset_y();

    data_out[i] = Pixel{x, y};
  }
}

//--------------------------------------------------------------------------
// Point

void Plot::add_point(boost::intrusive_ptr<Layer> const& layer,
    draw::PointStyle const& point_style,
    plot::Point const& plot_point)
{
  double x = plot_point.x();
  double y = plot_point.y();

  plot_point.create_draw_object({}, convert_x(x), convert_y(y), point_style);
  draw_layer_region_on(layer, plot_point.draw_object());
}

//--------------------------------------------------------------------------
// LinePiece

void Plot::apply_line_extend(double& x1, double& y1, double& x2, double& y2, LineExtend line_extend)
{
  if (line_extend != LineExtend::none)
  {
    double dx = x2 - x1;
    double dy = y2 - y1;
    double normal_x = dy;
    double normal_y = -dx;
    math::Hyperplane<2> line({normal_x, normal_y}, -(normal_x * x1 + normal_y * y1));
    math::Hyperblock<2> rectangle({range_[x_axis].min(), range_[y_axis].min()}, {range_[x_axis].max(), range_[y_axis].max()});
    auto intersections = rectangle.intersection_points(line);
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
}

void Plot::add_line(boost::intrusive_ptr<Layer> const& layer,
    draw::LineStyle const& line_style, LineExtend line_extend,
    LinePiece const& plot_line_piece)
{
  cairowindow::Point const& from = plot_line_piece.from();
  cairowindow::Point const& to = plot_line_piece.to();

  double x1 = from.x();
  double y1 = from.y();
  double x2 = to.x();
  double y2 = to.y();

  apply_line_extend(x1, y1, x2, y2, line_extend);
  plot_line_piece.draw_object_ = std::make_shared<draw::Line>(
      convert_x(x1), convert_y(y1), convert_x(x2), convert_y(y2),
      line_style);
  draw_layer_region_on(layer, plot_line_piece.draw_object_);
}

//--------------------------------------------------------------------------
// Connector

void Plot::add_connector(boost::intrusive_ptr<Layer> const& layer,
    draw::ConnectorStyle const& connector_style,
    Connector const& plot_connector)
{
  cairowindow::Point const& from = plot_connector.from();
  cairowindow::Point const& to = plot_connector.to();
  Connector::ArrowHeadShape arrow_head_shape_from = plot_connector.arrow_head_shape_from();
  Connector::ArrowHeadShape arrow_head_shape_to = plot_connector.arrow_head_shape_to();

  plot_connector.draw_object_ = std::make_shared<draw::Connector>(
      convert_x(from.x()), convert_y(from.y()), convert_x(to.x()), convert_y(to.y()),
      connector_style, arrow_head_shape_from, arrow_head_shape_to);
  if (need_print_)
    layer->start_printing_to(svg_cr_);
  layer->draw(plot_connector.draw_object_);
  plot_connector.draw_object_->draw_arrow_heads(layer);
  if (need_print_)
    layer->stop_printing();
}

//--------------------------------------------------------------------------
// Line

void Plot::add_line(boost::intrusive_ptr<Layer> const& layer,
    draw::LineStyle const& line_style,
    plot::Line const& plot_line)
{
  cairowindow::cs::Direction<csid::plot> const& direction = plot_line.direction();
  cairowindow::cs::Point<csid::plot> const& point = plot_line.point();

  double normal_x = -direction.y();
  double normal_y = direction.x();
  math::Hyperplane<2> line({normal_x, normal_y}, -(normal_x * point.x() + normal_y * point.y()));
  math::Hyperblock<2> rectangle({range_[x_axis].min(), range_[y_axis].min()}, {range_[x_axis].max(), range_[y_axis].max()});
  auto intersections = rectangle.intersection_points(line);

  // Is the line outside the plot area?
  if (intersections.empty())
    return;

  constexpr math::Hyperblock<2>::IntersectionPointIndex first{size_t{0}};
  constexpr math::Hyperblock<2>::IntersectionPointIndex second{size_t{1}};

  double x1 = intersections[first].coordinate(0);
  double y1 = intersections[first].coordinate(1);
  double x2 = intersections[second].coordinate(0);
  double y2 = intersections[second].coordinate(1);

  plot_line.create_draw_object({}, convert_x(x1), convert_y(y1), convert_x(x2), convert_y(y2), line_style);
  draw_layer_region_on(layer, plot_line.draw_object());
}

//--------------------------------------------------------------------------
// Rectangle

void Plot::add_rectangle(boost::intrusive_ptr<Layer> const& layer,
    draw::RectangleStyle const& rectangle_style, Rectangle const& plot_rectangle)
{
  double offset_x = plot_rectangle.offset_x();
  double offset_y = plot_rectangle.offset_y();
  double width = plot_rectangle.width();
  double height = plot_rectangle.height();

  plot_rectangle.draw_object_ = std::make_shared<draw::Rectangle>(
      convert_x(offset_x), convert_y(offset_y), convert_x(offset_x + width), convert_y(offset_y + height),
      rectangle_style);
  draw_layer_region_on(layer, plot_rectangle.draw_object_);
}

//--------------------------------------------------------------------------
// Circle

void Plot::add_circle(boost::intrusive_ptr<Layer> const& layer,
    draw::CircleStyle const& circle_style,
    Circle const& plot_circle)
{
  cairowindow::Point const& center = plot_circle.center();
  double radius = plot_circle.radius();

  plot_circle.draw_object_ = std::make_shared<draw::Circle>(
      cairowindow::Geometry{convert_x(center.x()), convert_y(center.y()), convert_x(radius) - convert_x(0), convert_y(0) - convert_y(radius)},
      circle_style);
  draw_layer_region_on(layer, plot_circle.draw_object_);
}

//--------------------------------------------------------------------------
// Arc

void Plot::add_arc(boost::intrusive_ptr<Layer> const& layer,
    draw::ArcStyle const& arc_style,
    Arc const& plot_arc)
{
  cairowindow::Point const& center = plot_arc.center();
  double radius = plot_arc.radius();
  double start_angle = plot_arc.start_angle();
  double end_angle = plot_arc.end_angle();

  plot_arc.draw_object_ = std::make_shared<draw::Arc>(
      convert_x(center.x()), convert_y(center.y()), -end_angle, -start_angle,
      std::max(convert_x(radius) - convert_x(0), convert_y(radius) - convert_y(0)),
      arc_style);
  draw_layer_region_on(layer, plot_arc.draw_object_);
}

//--------------------------------------------------------------------------
// BezierCurve

void Plot::add_bezier_curve(boost::intrusive_ptr<Layer> const& layer,
    draw::BezierCurveStyle const& bezier_curve_style,
    BezierCurve const& plot_bezier_curve)
{
  Vector const& P0 = plot_bezier_curve.P0();
  Vector const& C0 = plot_bezier_curve.C0();
  Vector const& C1 = plot_bezier_curve.C1();
  Vector const& P1 = plot_bezier_curve.P1();

  plot_bezier_curve.draw_object_ =
      std::make_shared<draw::BezierCurve>(
        convert_x(P0.x()), convert_y(P0.y()),
        convert_x(C0.x()), convert_y(C0.y()),
        convert_x(C1.x()), convert_y(C1.y()),
        convert_x(P1.x()), convert_y(P1.y()),
        bezier_curve_style);
  draw_layer_region_on(layer, plot_bezier_curve.draw_object_);
}

void Plot::add_bezier_curve_in_px(boost::intrusive_ptr<Layer> const& layer,
    draw::BezierCurveStyle const& bezier_curve_style,
    BezierCurve const& plot_bezier_curve_in_px)
{
  Vector const& P0 = plot_bezier_curve_in_px.P0();
  Vector const& C0 = plot_bezier_curve_in_px.C0();
  Vector const& C1 = plot_bezier_curve_in_px.C1();
  Vector const& P1 = plot_bezier_curve_in_px.P1();

  plot_bezier_curve_in_px.draw_object_ =
      std::make_shared<draw::BezierCurve>(
        P0.x(), P0.y(),
        C0.x(), C0.y(),
        C1.x(), C1.y(),
        P1.x(), P1.y(),
        bezier_curve_style);
  draw_layer_region_on(layer, plot_bezier_curve_in_px.draw_object_);
}

//--------------------------------------------------------------------------
// Curve

void Plot::add_bezier_fitter(boost::intrusive_ptr<Layer> const& layer,
    draw::LineStyle const& line_style,
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
    draw::PointStyle const& point_style,
#endif
    BezierFitter const& plot_bezier_fitter)
{
  std::vector<cairowindow::BezierCurve> const& bezier_curves = plot_bezier_fitter.result();

  // This call creates default constructed plot::BezierCurve objects from the "result" vector
  // of cairowindow::BezierCurve in bezier_curves.
  plot_bezier_fitter.draw_object_ = std::make_shared<draw::BezierFitter>(bezier_curves, line_style);
  for (BezierCurve const& plot_bezier_curve : plot_bezier_fitter.draw_object_->plot_bezier_curves())
    if (plot_bezier_curve.isfinite())
      add_bezier_curve(layer, line_style, plot_bezier_curve);
#if CAIROWINDOW_SHOW_BEZIER_CURVE_POINTS
  for (Point const& plot_bezier_curve_point : plot_bezier_fitter.draw_object_->plot_bezier_curve_points())
    add_point(layer, point_style, plot_bezier_curve_point);
#endif
}

//--------------------------------------------------------------------------
// Text

void Plot::add_text(boost::intrusive_ptr<Layer> const& layer,
    draw::TextStyle const& text_style,
    Text const& plot_text)
{
  cairowindow::Pixel position = plot_text.position();
  std::string const& text = plot_text.text();

  plot_text.draw_object_ = std::make_shared<draw::Text>(text, position.x(), position.y(), text_style);
  draw_layer_region_on(layer, plot_text.draw_object_);
}

//--------------------------------------------------------------------------
// Slider

Slider Plot::create_slider(boost::intrusive_ptr<Layer> const& layer,
    cairowindow::Geometry const& geometry, double start_value, double min_value, double max_value)
{
  Slider plot_slider(geometry, min_value, max_value);
  plot_slider.draw_object_ = std::make_shared<draw::Slider>(geometry.offset_x(), geometry.offset_y(), geometry.width(), geometry.height(),
      start_value, min_value, max_value);
  static_cast<draw::MultiRegion&>(*plot_slider.draw_object_).draw_regions_on(layer.get());
  Window* window = layer->window();
  std::shared_ptr<draw::Slider> slider_ptr = plot_slider.draw_object_;
  window->register_draggable(*this, slider_ptr.get(), [slider_ptr](cairowindow::Point const& point) -> cairowindow::Point { return point; });
  return plot_slider;
}

//--------------------------------------------------------------------------

void Plot::curve_to_bezier_curves(boost::intrusive_ptr<Layer> const& layer,
    Curve const& plot_curve, draw::BezierCurveStyle const& bezier_curve_style)
{
  std::vector<cairowindow::Point> const& plot_points = plot_curve.points();
  if (plot_points.size() < 2)
    return;
  std::vector<std::shared_ptr<draw::BezierCurve>>& bezier_curves = plot_curve.draw_object_->bezier_curves();
  if (plot_points.size() == 2)
  {
    // Draw a straight line: P(t) = P0 + t (P1 - P0)
    //   P(t) = B + V0 t + (A0/2) t² + J/6 t³
    Vector B{plot_points[0].raw()};
    Vector V0{plot_points[1].raw() - plot_points[0].raw()};
    BezierCurve bezier_curve(BezierCurveMatrix{{{B, V0, {}, {}}}});
    // Create and draw the draw object.
    add_bezier_curve(layer, bezier_curve_style, bezier_curve);
    // Move it into the Curve object.
    bezier_curves.emplace_back(std::move(bezier_curve.draw_object_));
    return;
  }
  if (plot_points.size() == 3)
  {
    // Implement.
    ASSERT(false);
  }

//  cairowindow::Rectangle const& g = plot_area_.geometry();
  cairowindow::Point prev_point;

  for (int i = 1; i < plot_points.size(); ++i)
  {
    Vector B{plot_points[i - 1].raw()};
    Vector V0{plot_points[i].raw() - plot_points[i - 1].raw()};
    BezierCurve bezier_curve(BezierCurveMatrix{{{B, V0, {}, {}}}});
    add_bezier_curve(layer, bezier_curve_style, bezier_curve);
    bezier_curves.emplace_back(std::move(bezier_curve.draw_object_));
  }
}

Curve Plot::create_curve(boost::intrusive_ptr<Layer> const& layer,
  draw::BezierCurveStyle const& bezier_curve_style,
  std::vector<cairowindow::Point>&& points)
{
  Curve plot_curve(std::move(points), std::make_shared<draw::Curve>(bezier_curve_style));
  curve_to_bezier_curves(layer, plot_curve, bezier_curve_style);
  return plot_curve;
}

cairowindow::Geometry Plot::update_grabbed(utils::Badge<Window>, ClickableIndex grabbed_point, double pixel_x, double pixel_y)
{
  Draggable* draggable = draggables_[grabbed_point];
  // If convert is not true then pixel_x, pixel_y are actually cairowindow::Point coordinates (aka csid::plot).
  cairowindow::Point new_position = draggable->convert() ? convert_from_pixel(Pixel{pixel_x, pixel_y}) : cairowindow::Point{pixel_x, pixel_y};

  if (draggable_restrictions_[grabbed_point])
    new_position = draggable_restrictions_[grabbed_point](new_position);

  draggable->moved(this, new_position);

  return draggable->geometry();
}

void Plot::apply_restrictions(utils::Badge<Window>, ClickableIndex clickable_index, cairowindow::Point& new_position)
{
  if (draggable_restrictions_[clickable_index])
    new_position = draggable_restrictions_[clickable_index](new_position);
}

void Slider::set_value(double value)
{
  min_value_ = std::min(min_value_, value);
  max_value_ = std::max(max_value_, value);
  draw_object_->set_value(value, min_value_, max_value_);
}

} // namespace cairowindow::plot
