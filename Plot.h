#pragma once

#include "draw/PlotArea.h"
#include "draw/Text.h"
#include "Range.h"
#include <boost/intrusive_ptr.hpp>
#include <string>

namespace cairowindow::plot {

class Plot
{
 private:
  draw::PlotArea plot_area_;
  draw::Text title_;
  draw::Text xlabel_;
  draw::Text ylabel_;
  Range x_range_;
  Range y_range_;

  struct TitleStyleDefaults : draw::DefaultTextStyleDefaults
  {
    static constexpr draw::TextPosition position = draw::centered;
    static constexpr double font_size = 24.0;
    static constexpr double offset = 0.0;
  };

  struct XLabelStyleDefaults : draw::DefaultTextStyleDefaults
  {
    static constexpr draw::TextPosition position = draw::centered_below;
    static constexpr double font_size = 18.0;
  };

  struct YLabelStyleDefaults : draw::DefaultTextStyleDefaults
  {
    static constexpr draw::TextPosition position = draw::centered_above;
    static constexpr double font_size = 18;
    static constexpr double rotation = -0.5 * M_PI;
  };

 public:
  using TitleStyle = draw::TextStyle<TitleStyleDefaults>;
  using XLabelStyle = draw::TextStyle<XLabelStyleDefaults>;
  using YLabelStyle = draw::TextStyle<YLabelStyleDefaults>;

 public:
  Plot(Rectangle const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, TitleStyle title_style) :
    plot_area_(axes_geometry(geometry), plot_area_style),
    title_(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - title_style.offset, title_style) { }

  Plot(Rectangle const& geometry, draw::PlotAreaStyle plot_area_style, std::string title, TitleStyle title_style,
      std::string xlabel, XLabelStyle xlabel_style, std::string ylabel, YLabelStyle ylabel_style) :
    plot_area_(axes_geometry(geometry), plot_area_style),
    title_(title, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() - 0.5 * plot_area_.geometry().offset_y() - title_style.offset, title_style),
    xlabel_(xlabel, plot_area_.geometry().offset_x() + 0.5 * plot_area_.geometry().width(),
        plot_area_.geometry().offset_y() + plot_area_.geometry().height() + xlabel_style.offset, xlabel_style),
    ylabel_(ylabel, plot_area_.geometry().offset_x() - ylabel_style.offset,
        plot_area_.geometry().offset_y() + 0.5 * plot_area_.geometry().height(), ylabel_style) { }

  Plot(Rectangle const& geometry, draw::PlotAreaStyle plot_area_style) :
    plot_area_(axes_geometry(geometry), plot_area_style) { }

  void set_xrange(Range x_range)
  {
    x_range_ = x_range;
  }

  void set_yrange(Range y_range)
  {
    y_range_ = y_range;
  }

  void add_to(boost::intrusive_ptr<Layer> const& layer);

 private:
  Rectangle axes_geometry(Rectangle const& geometry);
};

} // namespace cairowindow::plot
