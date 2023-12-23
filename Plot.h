#pragma once

#include "draw/PlotArea.h"
#include "draw/Text.h"
#include "Range.h"
#include <boost/intrusive_ptr.hpp>
#include <string>
#include <vector>

namespace cairowindow::plot {

class Plot
{
  static constexpr int number_of_axes = draw::PlotArea::number_of_axes;

 private:
  draw::PlotArea plot_area_;
  draw::Text title_;
  draw::Text xlabel_;
  draw::Text ylabel_;
  std::array<Range, number_of_axes> range_;
  std::array<int, number_of_axes> range_ticks_{{10, 10}};
  std::array<std::vector<std::unique_ptr<draw::Text>>, number_of_axes> labels;

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
    static constexpr double offset = -10.0;
  };

  struct YLabelStyleDefaults : LabelStyleDefaults
  {
    static constexpr draw::TextPosition position = draw::centered_above;
    static constexpr double rotation = -0.5 * M_PI;
    static constexpr double offset = -XLabelStyleDefaults::offset;
  };

 public:
  using TitleStyle = draw::TextStyle<TitleStyleDefaults>;
  using XLabelStyle = draw::TextStyle<XLabelStyleDefaults>;
  using YLabelStyle = draw::TextStyle<YLabelStyleDefaults>;
  using LabelStyle = draw::TextStyle<LabelStyleDefaults>;

 public:
  Plot(Rectangle const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, TitleStyle title_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width), plot_area_style),
    title_(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - title_style.offset, title_style) { }

  Plot(Rectangle const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, TitleStyle title_style,
      std::string xlabel, XLabelStyle xlabel_style, std::string ylabel, YLabelStyle ylabel_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width), plot_area_style),
    title_(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_style.offset, title_style),
    xlabel_(xlabel, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() + plot_area_.geometry().height() - XLabelStyleDefaults::offset, xlabel_style),
    ylabel_(ylabel, plot_area_.geometry().offset_x() - YLabelStyleDefaults::offset,
        plot_area_.geometry().offset_y() + 0.5 * plot_area_.geometry().height(), ylabel_style) { }

  Plot(Rectangle const& geometry, draw::PlotAreaStyle plot_area_style) :
    plot_area_(axes_geometry(geometry, plot_area_style.axes_line_width), plot_area_style) { }

  void set_range(int axis, Range range)
  {
    DoutEntering(dc::notice, "Plot::set_range(" << axis << ", " << range << ")");
    range_[axis] = range;
    range_ticks_[axis] = draw::PlotArea::calculate_range_ticks(range_[axis]);
    Dout(dc::notice, "range_[" << axis << "] = " << range_[axis] << "; range_ticks_[" << axis << "] = " << range_ticks_[axis]);
  }

  void set_xrange(Range x_range) { set_range(x_axis, x_range); }
  void set_yrange(Range y_range) { set_range(y_axis, y_range); }

  void add_to(boost::intrusive_ptr<Layer> const& layer);

 private:
  Rectangle axes_geometry(Rectangle const& geometry, double axes_line_width);
};

} // namespace cairowindow::plot
