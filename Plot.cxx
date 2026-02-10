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

void Plot::update_plot_transform_pixels()
{
  if (range_[x_axis].size() == 0.0 || range_[y_axis].size() == 0.0)
    return;

  // Initialize cs_transform_pixels_.
  // x' = (x - xmin) / xrange * width + offset_x
  // y' = (ymax - y) / yrange * height + offset_y
  cairowindow::Geometry const& g = plot_area_.geometry();
  double const sx = g.width() / range_[x_axis].size();
  double const sy = g.height() / range_[y_axis].size();
  double const tx = g.offset_x() - range_[x_axis].min() * sx;
  double const ty = g.offset_y() + range_[y_axis].max() * sy;
  cs_transform_pixels_ = math::Transform<csid::plot, csid::pixels>{}
    .translate(math::TranslationVector<csid::pixels>::create_from_cs_values(tx, ty))
    .scale(sx, -sy);
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

  update_plot_transform_pixels();

  // Register this plot and its geometry with the associated Window so that we can find which printable is under the mouse if needed.
  layer->window()->add_printable(this);
}

//--------------------------------------------------------------------------
// LinePiece

void Plot::add_line_piece(boost::intrusive_ptr<Layer> const& layer,
    draw::LineStyle const& line_style, LineExtend line_extend,
    LinePiece const& plot_line_piece)
{
  math::Hyperblock<2> clip_rectangle({range_[x_axis].min(), range_[y_axis].min()}, {range_[x_axis].max(), range_[y_axis].max()});
  add_clipped_line_piece(layer, line_style, line_extend, plot_line_piece, clip_rectangle);
}

//--------------------------------------------------------------------------
// Line

void Plot::add_line(boost::intrusive_ptr<Layer> const& layer,
    draw::LineStyle const& line_style,
    plot::Line const& plot_line)
{
  math::Hyperblock<2> clip_rectangle({range_[x_axis].min(), range_[y_axis].min()}, {range_[x_axis].max(), range_[y_axis].max()});
  add_clipped_line(layer, line_style, plot_line, clip_rectangle);
}

//--------------------------------------------------------------------------
// Circle

//--------------------------------------------------------------------------
// BezierCurve

void Plot::add_bezier_curve(boost::intrusive_ptr<Layer> const& layer,
    draw::BezierCurveStyle const& bezier_curve_style,
    BezierCurve const& plot_bezier_curve)
{
  // Use explicit conversion from math::Point<2> to cairowindow::Point (aka math::cs::Point<csid::plot>).
  cairowindow::Point const P0{plot_bezier_curve.P0()};
  cairowindow::Point const C0{plot_bezier_curve.C0()};
  cairowindow::Point const C1{plot_bezier_curve.C1()};
  cairowindow::Point const P1{plot_bezier_curve.P1()};

  math::cs::Point<csid::pixels> const P0_pixels = P0 * cs_transform_pixels_;
  math::cs::Point<csid::pixels> const C0_pixels = C0 * cs_transform_pixels_;
  math::cs::Point<csid::pixels> const C1_pixels = C1 * cs_transform_pixels_;
  math::cs::Point<csid::pixels> const P1_pixels = P1 * cs_transform_pixels_;

  plot_bezier_curve.draw_object_ =
      std::make_shared<draw::BezierCurve>(
        P0_pixels.x(), P0_pixels.y(),
        C0_pixels.x(), C0_pixels.y(),
        C1_pixels.x(), C1_pixels.y(),
        P1_pixels.x(), P1_pixels.y(),
        bezier_curve_style);
  draw_layer_region_on(layer, plot_bezier_curve.draw_object_);
}

void Plot::add_bezier_curve_in_px(boost::intrusive_ptr<Layer> const& layer,
    draw::BezierCurveStyle const& bezier_curve_style,
    BezierCurve const& plot_bezier_curve_in_px)
{
  // This BezierCurve is in pixels coordinates.
  math::Point<2> const& P0_pixels = plot_bezier_curve_in_px.P0();
  math::Point<2> const& C0_pixels = plot_bezier_curve_in_px.C0();
  math::Point<2> const& C1_pixels = plot_bezier_curve_in_px.C1();
  math::Point<2> const& P1_pixels = plot_bezier_curve_in_px.P1();

  plot_bezier_curve_in_px.draw_object_ =
      std::make_shared<draw::BezierCurve>(
        P0_pixels.x(), P0_pixels.y(),
        C0_pixels.x(), C0_pixels.y(),
        C1_pixels.x(), C1_pixels.y(),
        P1_pixels.x(), P1_pixels.y(),
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
  std::shared_ptr<draw::Slider> slider_ptr = plot_slider.draw_object_;  // Kept alive by the 'restriction' lambda below.
  window->register_draggable(slider_ptr.get(),
      [slider_ptr](math::cs::Point<csid::pixels> const& position_pixels) -> math::cs::Point<csid::pixels>
      {
        return position_pixels;
      });
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
    math::Vector<2> B{plot_points[0].raw()};
    math::Vector<2> V0{plot_points[1].raw() - plot_points[0].raw()};
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
    math::Vector<2> B{plot_points[i - 1].raw()};
    math::Vector<2> V0{plot_points[i].raw() - plot_points[i - 1].raw()};
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

void Slider::set_value(double value)
{
  min_value_ = std::min(min_value_, value);
  max_value_ = std::max(max_value_, value);
  draw_object_->set_value(value, min_value_, max_value_);
}

} // namespace cairowindow::plot
