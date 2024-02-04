#include "sys.h"
#include "Layer.h"
#include "Plot.h"
#include "intersection_points.h"
#include "utils/square.h"
#include <limits>
#include <cmath>
#include <iomanip>

namespace cairowindow::plot {

// Calculate the axes geometry from the full plot geometry.
cairowindow::Rectangle Plot::axes_geometry(cairowindow::Rectangle const& geometry, double axes_line_width)
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
    cairowindow::Rectangle const& geometry = plot_area_.geometry();
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
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_->style().offset);
    layer->draw(title_);
  }

  // Set ranges on the plot area and draw it.
  for (int axis = x_axis; axis <= y_axis; ++axis)
    plot_area_.set_range(axis, range_[axis].min(), range_[axis].max());
  layer->draw(&plot_area_);

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
            x, y, LabelStyle{.position = (axis == x_axis) ? draw::centered_below : draw::centered_left_of}));
      layer->draw(labels_[axis].back());
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
    layer->draw(xlabel_);
  }

  if (ylabel_)
  {
    ylabel_->rel_move_to(-max_width, 0);
    layer->draw(ylabel_);
  }
}

double Plot::convert_x(double x) const
{
  cairowindow::Rectangle const& g = plot_area_.geometry();
  x = (x - range_[x_axis].min()) / range_[x_axis].size();
  x *= g.width();
  x += g.offset_x();
  return x;
}

double Plot::convert_y(double y) const
{
  cairowindow::Rectangle const& g = plot_area_.geometry();
  y = (range_[y_axis].max() - y) / range_[y_axis].size();
  y *= g.height();
  y += g.offset_y();
  return y;
}

Pixel Plot::convert_to_pixel(cairowindow::Point const& point) const
{
  double x = point.x();
  double y = point.y();

  cairowindow::Rectangle const& g = plot_area_.geometry();

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

//--------------------------------------------------------------------------
// Point

void Plot::add_point(boost::intrusive_ptr<Layer> const& layer,
    draw::PointStyle const& point_style,
    Point const& plot_point)
{
  double x = plot_point.x();
  double y = plot_point.y();

  plot_point.draw_object_ = std::make_shared<draw::Point>(
      convert_x(x), convert_y(y),
      point_style);
  layer->draw(plot_point.draw_object_);
}

//--------------------------------------------------------------------------
// LinePiece

void Plot::apply_line_extend(double& x1, double& y1, double& x2, double& y2, LineExtend line_extend)
{
  if (line_extend != LineExtend::none)
  {
    double dx = x2 - x1;
    double dy = y2 - y1;
    intersections::HyperPlane<double, 2> line({dy, -dx}, y1 * dx - x1 * dy);
    intersections::HyperBlock<double, 2> rectangle({range_[x_axis].min(), range_[y_axis].min()}, {range_[x_axis].max(), range_[y_axis].max()});
    auto intersections = rectangle.intersection_points(line);
    if (!intersections.empty())
    {
      if (line_extend != LineExtend::to)
      {
        x1 = intersections[0][0];
        y1 = intersections[0][1];
      }
      if (line_extend != LineExtend::from)
      {
        x2 = intersections[1][0];
        y2 = intersections[1][1];
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
  layer->draw(plot_line_piece.draw_object_);
}

//--------------------------------------------------------------------------
// Connector

void Plot::add_connector(boost::intrusive_ptr<Layer> const& layer,
    draw::LineStyle const& line_style, Color fill_color,
    Connector::ArrowHeadShape arrow_head_shape_from, Connector::ArrowHeadShape arrow_head_shape_to,
    Connector const& plot_connector)
{
  cairowindow::Point const& from = plot_connector.from();
  cairowindow::Point const& to = plot_connector.to();

  plot_connector.draw_object_ = std::make_shared<draw::Connector>(
      convert_x(from.x()), convert_y(from.y()), convert_x(to.x()), convert_y(to.y()),
      line_style, fill_color, arrow_head_shape_from, arrow_head_shape_to);
  layer->draw(plot_connector.draw_object_);
  plot_connector.draw_object_->draw_arrow_heads(layer);
}

//--------------------------------------------------------------------------
// Line

void Plot::add_line(boost::intrusive_ptr<Layer> const& layer,
    draw::LineStyle const& line_style,
    Line const& plot_line)
{
  cairowindow::Direction const& direction = plot_line.direction();
  cairowindow::Point const& point = plot_line.point();

  double normal_x = -direction.y();
  double normal_y = direction.x();
  intersections::HyperPlane<double, 2> line({normal_x, normal_y}, -point.x() * normal_x - point.y() * normal_y);
  intersections::HyperBlock<double, 2> rectangle({range_[x_axis].min(), range_[y_axis].min()}, {range_[x_axis].max(), range_[y_axis].max()});
  auto intersections = rectangle.intersection_points(line);

  // Is the line outside the plot area?
  if (intersections.empty())
    return;

  double x1 = intersections[0][0];
  double y1 = intersections[0][1];
  double x2 = intersections[1][0];
  double y2 = intersections[1][1];

  plot_line.draw_object_ = std::make_shared<draw::Line>(
      convert_x(x1), convert_y(y1), convert_x(x2), convert_y(y2),
      line_style);
  layer->draw(plot_line.draw_object_);
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
  layer->draw(plot_rectangle.draw_object_);
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
      cairowindow::Rectangle{convert_x(center.x()), convert_y(center.y()), convert_x(radius) - convert_x(0), convert_y(0) - convert_y(radius)},
      circle_style);
  layer->draw(plot_circle.draw_object_);
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
  layer->draw(plot_arc.draw_object_);
}

//--------------------------------------------------------------------------
// BezierCurve

void Plot::add_bezier_curve(boost::intrusive_ptr<Layer> const& layer,
    draw::BezierCurveStyle const& bezier_curve_style,
    BezierCurve const& plot_bezier_curve)
{
  Vector const& P0 = plot_bezier_curve.P0();
  Vector const& C1 = plot_bezier_curve.C1();
  Vector const& C2 = plot_bezier_curve.C2();
  Vector const& P1 = plot_bezier_curve.P1();

  plot_bezier_curve.draw_object_ =
      std::make_shared<draw::BezierCurve>(
        convert_x(P0.x()), convert_y(P0.y()),
        convert_x(C1.x()), convert_y(C1.y()),
        convert_x(C2.x()), convert_y(C2.y()),
        convert_x(P1.x()), convert_y(P1.y()),
        bezier_curve_style);
  layer->draw(plot_bezier_curve.draw_object_);
}

//--------------------------------------------------------------------------
// Text

void Plot::add_text(boost::intrusive_ptr<Layer> const& layer,
    draw::TextStyle<> const& text_style,
    Text const& plot_text)
{
  cairowindow::Pixel position = plot_text.position();
  std::string const& text = plot_text.text();

  plot_text.draw_object_ = std::make_shared<draw::Text>(text, position.x(), position.y(), text_style);
  layer->draw(plot_text.draw_object_);
}

//--------------------------------------------------------------------------
// Slider

Slider Plot::create_slider(boost::intrusive_ptr<Layer> const& layer,
    cairowindow::Rectangle const& geometry, double start_value, double min_value, double max_value)
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

void Plot::curve_to_lines(boost::intrusive_ptr<Layer> const& layer, Curve const& plot_curve, draw::LineStyle const& line_style)
{
  cairowindow::Rectangle const& g = plot_area_.geometry();
  double prev_point_x = std::numeric_limits<double>::quiet_NaN();
  double prev_point_y = std::numeric_limits<double>::quiet_NaN();
  int count = 0;
  std::vector<std::shared_ptr<draw::Line>>& lines = plot_curve.draw_object_->lines();
  for (cairowindow::Point const& point : plot_curve.points())
  {
    double x = convert_x(point.x());
    double y = convert_y(point.y());
    if (!std::isnan(prev_point_x) && !std::isnan(prev_point_y) && !std::isnan(x) && !std::isnan(y))
    {
      lines.emplace_back(std::make_shared<draw::Line>(prev_point_x, prev_point_y, x, y, line_style));
      layer->draw(lines.back());
    }
    prev_point_x = x;
    prev_point_y = y;
    ++count;
  }
}

Curve Plot::create_curve(boost::intrusive_ptr<Layer> const& layer,
    draw::LineStyle const& line_style,
    std::vector<cairowindow::Point>&& points)
{
  Curve plot_curve(std::move(points), std::make_shared<draw::Curve>(line_style));
  curve_to_lines(layer, plot_curve, line_style);
  return plot_curve;
}

cairowindow::Rectangle Plot::update_grabbed(utils::Badge<Window>, ClickableIndex grabbed_point, double pixel_x, double pixel_y)
{
  double x = pixel_x;
  double y = pixel_y;
  Draggable* draggable = draggables_[grabbed_point];

  if (draggable->convert())
  {
    cairowindow::Rectangle const& g = plot_area_.geometry();
    x -= g.offset_x();
    x /= g.width();
    x *= range_[x_axis].size();
    x += range_[x_axis].min();

    y -= g.offset_y();
    y /= g.height();
    y *= range_[y_axis].size();
    y = range_[y_axis].max() - y;
  }

  cairowindow::Point new_position{x, y};

  if (draggable_restrictions_[grabbed_point])
    new_position = draggable_restrictions_[grabbed_point](new_position);

  draggable->moved(this, new_position);

  return draggable->geometry();
}

void Point::moved(Plot* plot, cairowindow::Point const& new_position)
{
  Point plot_point(new_position);
  plot->add_point(draw_object_->layer(), draw_object_->point_style(), plot_point);
  *this = plot_point;
}

void Slider::set_value(double value)
{
  min_value_ = std::min(min_value_, value);
  max_value_ = std::max(max_value_, value);
  draw_object_->set_value(value, min_value_, max_value_);
}

} // namespace cairowindow::plot
