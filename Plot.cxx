#include "sys.h"
#include "Layer.h"
#include "Plot.h"
#include "intersection_points.h"
#include <cmath>
#include <iomanip>

namespace cairowindow::plot {

// Calculate the axes geometry from the full plot geometry.
Rectangle Plot::axes_geometry(Rectangle const& geometry, double axes_line_width)
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
    Rectangle const& geometry = plot_area_.geometry();
    double pixels_per_x_unit = geometry.width() / range_[x_axis].size();
    double pixels_per_y_unit = geometry.height() / range_[y_axis].size();
    double required_scale = pixels_per_y_unit / pixels_per_x_unit;
    if (required_scale < 1.0)
      plot_area_.set_geometry({ geometry.offset_x(), geometry.offset_y(), geometry.width() * required_scale, geometry.height() });
    else
      plot_area_.set_geometry({ geometry.offset_x(), geometry.offset_y(), geometry.width(), geometry.height() / required_scale });
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
  Rectangle const& g = plot_area_.geometry();
  x = (x - range_[x_axis].min()) / range_[x_axis].size();
  x *= g.width();
  x += g.offset_x();
  return x;
}

double Plot::convert_y(double y) const
{
  Rectangle const& g = plot_area_.geometry();
  y = (range_[y_axis].max() - y) / range_[y_axis].size();
  y *= g.height();
  y += g.offset_y();
  return y;
}

void Plot::add_point(boost::intrusive_ptr<Layer> const& layer,
    Point const& plot_point, draw::PointStyle const& point_style)
{
  plot_point.draw_object_ = std::make_shared<draw::Point>(convert_x(plot_point.x()), convert_y(plot_point.y()), point_style);
  layer->draw(plot_point.draw_object_);
}

void Plot::add_text(boost::intrusive_ptr<Layer> const& layer,
    Text const& text, draw::TextStyle<> const& text_style)
{
  text.draw_object_ = std::make_shared<draw::Text>(text.text(), convert_x(text.position().x()), convert_y(text.position().y()), text_style);
  layer->draw(text.draw_object_);
}

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
    LinePiece const& plot_line_piece, draw::LineStyle const& line_style, LineExtend line_extend)
{
  double x1 = plot_line_piece.from().x();
  double y1 = plot_line_piece.from().y();
  double x2 = plot_line_piece.to().x();
  double y2 = plot_line_piece.to().y();

  apply_line_extend(x1, y1, x2, y2, line_extend);
  plot_line_piece.draw_object_ = std::make_shared<draw::Line>(convert_x(x1), convert_y(y1), convert_x(x2), convert_y(y2), line_style);
  layer->draw(plot_line_piece.draw_object_);
}

void Plot::add_connector(boost::intrusive_ptr<Layer> const& layer,
    Connector const& plot_connector, draw::LineStyle const& line_style, Color fill_color)
{
  double x1 = plot_connector.from().x();
  double y1 = plot_connector.from().y();
  double x2 = plot_connector.to().x();
  double y2 = plot_connector.to().y();

  plot_connector.draw_object_ = std::make_shared<draw::Connector>(convert_x(x1), convert_y(y1), convert_x(x2), convert_y(y2), line_style,
      fill_color, plot_connector.arrow_head_shape_from(), plot_connector.arrow_head_shape_to());
  layer->draw(plot_connector.draw_object_);
  plot_connector.draw_object_->draw_arrow_heads(layer);
}

Line Plot::create_line(boost::intrusive_ptr<Layer> const& layer,
    cairowindow::Point const& point, cairowindow::Direction const& direction, draw::LineStyle const& line_style)
{
  double normal_x = -direction.y();
  double normal_y = direction.x();
  intersections::HyperPlane<double, 2> line({normal_x, normal_y}, -point.x() * normal_x - point.y() * normal_y);
  intersections::HyperBlock<double, 2> rectangle({range_[x_axis].min(), range_[y_axis].min()}, {range_[x_axis].max(), range_[y_axis].max()});
  auto intersections = rectangle.intersection_points(line);

  // Is the line outside the plot area?
  if (intersections.empty())
    return {direction, point};

  double x1 = intersections[0][0];
  double y1 = intersections[0][1];
  double x2 = intersections[1][0];
  double y2 = intersections[1][1];

  Line plot_line(direction, point, std::make_shared<draw::Line>(convert_x(x1), convert_y(y1), convert_x(x2), convert_y(y2), line_style));
  layer->draw(plot_line.draw_object_);
  return plot_line;
}

Point Plot::create_point(boost::intrusive_ptr<Layer> const& layer,
    cairowindow::Point const& point, draw::PointStyle const& point_style)
{
  Point plot_point(point, std::make_shared<draw::Point>(convert_x(point.x()), convert_y(point.y()), point_style));
  layer->draw(plot_point.draw_object_);
  return plot_point;
}

Circle Plot::create_circle(boost::intrusive_ptr<Layer> const& layer,
    cairowindow::Point const& center, double radius, draw::CircleStyle const& circle_style)
{
  Circle plot_circle(center, radius, std::make_shared<draw::Circle>(
        Rectangle{convert_x(center.x()), convert_y(center.y()), convert_x(radius) - convert_x(0), convert_y(0) - convert_y(radius)},
        circle_style));
  layer->draw(plot_circle.draw_object_);
  return plot_circle;
}

Text Plot::create_text(boost::intrusive_ptr<Layer> const& layer,
    cairowindow::Point const& position, std::string const& text, draw::TextStyle<> const& text_style)
{
  Text plot_text(position, text, std::make_shared<draw::Text>(text, convert_x(position.x()), convert_y(position.y()), text_style));
  layer->draw(plot_text.draw_object_);
  return plot_text;
}

LinePiece Plot::create_line(boost::intrusive_ptr<Layer> const& layer,
    cairowindow::Point const& from, cairowindow::Point const& to, draw::LineStyle const& line_style, LineExtend line_extend)
{
  double x1 = from.x();
  double y1 = from.y();
  double x2 = to.x();
  double y2 = to.y();

  apply_line_extend(x1, y1, x2, y2, line_extend);
  LinePiece plot_line_piece(from, to, std::make_shared<draw::Line>(convert_x(x1), convert_y(y1), convert_x(x2), convert_y(y2), line_style));
  layer->draw(plot_line_piece.draw_object_);
  return plot_line_piece;
}

Curve Plot::create_curve(boost::intrusive_ptr<Layer> const& layer,
    std::vector<cairowindow::Point>&& points, draw::LineStyle const& line_style)
{
  Curve plot_curve(std::move(points), std::make_shared<draw::Curve>(line_style));

  Rectangle const& g = plot_area_.geometry();
  double prev_point_x;
  double prev_point_y;
  int count = 0;
  std::vector<std::shared_ptr<draw::Line>>& lines = plot_curve.draw_object_->lines();
  for (cairowindow::Point const& point : plot_curve.points())
  {
    double x = convert_x(point.x());
    double y = convert_y(point.y());
    if (count > 0)
    {
      lines.emplace_back(std::make_shared<draw::Line>(prev_point_x, prev_point_y, x, y, line_style));
      layer->draw(lines.back());
    }
    prev_point_x = x;
    prev_point_y = y;
    ++count;
  }

  return plot_curve;
}

} // namespace cairowindow::plot
