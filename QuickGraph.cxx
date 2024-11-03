#include "sys.h"
#include "QuickGraph.h"
#include "BezierFitter.h"
#include "debug.h"

namespace cairowindow {

QuickGraph::QuickGraph(std::string const& title, std::string const& x_label, std::string const& y_label, Range x_range) :
  title_(title), x_label_(x_label), y_label_(y_label), x_range_(x_range),
  window_(title, 1200, 900),
  background_layer_(window_.create_background_layer<Layer>(color::white COMMA_DEBUG_ONLY("QuickGraph::background_layer_"))),
  second_layer_(window_.create_layer<Layer>({} COMMA_DEBUG_ONLY("QuickGraph::second_layer_"))),
  event_loop_([&](){
    Debug(NAMESPACE_DEBUG::init_thread("QuickGraph::event_loop_"));
    EventLoop event_loop = window_.run();
    event_loop.set_cleanly_terminated();
  }),
  plot_(window_.geometry(), { .grid = {.color = color::orange} },
      title, {},
      x_label, {},
      y_label, {})
{
  DoutEntering(dc::notice, "QuickGraph::QuickGraph(\"" <<
      title << "\", \"" << x_label << "\", \"" << y_label << "\", " << x_range << ") [" << this << "]");
}

void QuickGraph::add_function(std::function<double(double)> const& f, draw::LineStyle const& line_style)
{
  DoutEntering(dc::notice, "QuickGraph::add_function(f, " << line_style << ") [" << this << "]");

  Rectangle viewport;
  if (empty_)
  {
    // Make a very rough estimate of the y-axis range by evaluating the function f at a few points.
    double const width = x_range_.size();
    double y_min, y_max;
    y_min = y_max = f(x_range_.min());
    for (int i = 1; i <= 8; ++i)
    {
      double x_probe = i * x_range_.size() / 8 + x_range_.min();
      double y_probe = f(x_probe);
      if (y_probe < y_min)
        y_min = y_probe;
      if (y_probe > y_max)
        y_max = y_probe;
    }

    // Construct a reasonable viewport from this (things outside the viewport are not drawn).
    double probed_height = y_max - y_min;
    y_min -= 0.5 * probed_height;
    viewport = Rectangle{x_range_.min(), y_min, width, 2 * probed_height};
  }
  else
    viewport = plot_.viewport();

  // Fit piecewise Bezier curves to the function.
  plot_bezier_fitter_.emplace_back([&f](double x) -> Point { return {x, f(x)}; }, x_range_, viewport);

  if (empty_)
  {
    // Get a much more accurate estimate of the y-axis range from the fitted Bezier curves.
    Range y_range;
    std::vector<BezierCurve> const& bezier_curves = plot_bezier_fitter_.back().result();
    bool first = true;
    for (BezierCurve const& bezier_curve : bezier_curves)
    {
      Rectangle extents = bezier_curve.extents();
      if (first)
        y_range = Range{extents.offset_y(), extents.offset_y() + extents.height()};
      else
        y_range = Range{std::min(y_range.min(), extents.offset_y()), std::max(y_range.max(), extents.offset_y() + extents.height())};
      first = false;
    }
    Dout(dc::notice, "y_range = " << y_range);

    plot_.set_xrange(x_range_);
    plot_.set_yrange(y_range);
    plot_.add_to(background_layer_, false);

    empty_ = false;
  }

  // Suppress immediate updating of the window for each created item, in order to avoid flickering.
  window_.set_send_expose_events(false);

  // Draw bezier_fitter.
  plot_.add_bezier_fitter(second_layer_, line_style, plot_bezier_fitter_.back());

  // Flush all expose events related to the drawing done above.
  window_.set_send_expose_events(true);
}

void QuickGraph::add_point(Point P, draw::PointStyle const& point_style)
{
  plot_points_.emplace_back(P);
  plot_.add_point(second_layer_, point_style, plot_points_.back());
}

void QuickGraph::wait_for_keypress()
{
  // Block until a key press.
  window_.handle_input_events();
  event_loop_.join();
}

QuickGraph::~QuickGraph()
{
  if (!window_.destroyed())
  {
    window_.close();
    event_loop_.join();
  }
}

} // namespace cairowindow
