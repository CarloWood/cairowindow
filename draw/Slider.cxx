#include "sys.h"
#include "Slider.h"
#include "cairowindow/Layer.h"
#include <string>
#include <sstream>
#include "debug.h"

namespace cairowindow::draw {

SliderHandle::SliderHandle(Slider const* parent) :
  Shape(
     {parent->track().offset_x() + ((parent->track().orientation() == SliderTrack::horizontal) ? parent->horizontal_position_pixels() : 0),
      parent->track().offset_y() + ((parent->track().orientation() == SliderTrack::vertical) ? parent->vertical_position_pixels() : 0),
      1.5 * parent->track().width(), 1.5 * parent->track().width()},
     { .line_color = Color{0.5, 0.5, 0.5, 0.5},
       .fill_color = color::white,
       .line_width = SliderTrack::track_width,
       .position = at_corner,
       .shape = ellipse })
{
  DoutEntering(dc::cairowindow, "SliderHandle(" << parent << ") [" << this << "]");
}

void SliderTrack::draw_highlighted_track(cairo_t* cr, Slider const* parent) const
{
#ifdef CWDEBUG
  using namespace debugcairo;
#endif
  if (orientation_ == horizontal)
    cairo_rectangle(cr, offset_x_ + stop_offset, offset_y_ - 0.5 * track_width, parent->horizontal_position_pixels() - stop_offset, track_width);
  else
    cairo_rectangle(cr, offset_x_ - 0.5 * track_width, offset_y_ + parent->vertical_position_pixels(),
        track_width, length_ - stop_offset - parent->vertical_position_pixels());
}

struct DrawHighlightedTrack
{
  Slider* slider;
  StrokeExtents operator()(cairo_t* cr)
  {
    slider->do_set_fill_color(cr);
    slider->track().draw_highlighted_track(cr, slider);
    return slider->do_fill(cr);
  }
};

void Slider::set_value_text()
{
  std::ostringstream value_str;
  value_str << std::setprecision(3) << (rel_value_ * (max_value_ - min_value_) + min_value_);
  //FIXME: position correctly for horizontal sliders.
  value_text_ = std::make_shared<Text>(value_str.str(), track_->offset_x(), track_->offset_y(),
      TextStyle<>{.position = centered_above, .font_size = 14, .color = color::black, .font_family = "sans-serif", .offset = 10, .rotation = 0});
}

void Slider::draw_regions_on(Layer* layer)
{
  DoutEntering(dc::cairowindow, "Slider::draw_regions_on(" << layer << ")");

  highlighted_track_ = std::make_shared<LayerRegion>(DrawHighlightedTrack{this});
  set_value_text();

  layer->draw(track_);
  layer->draw(value_text_);
  layer->draw(highlighted_track_);
  layer->draw(handle_);
}

void Slider::moved(plot::Plot* plot, cairowindow::Point const& new_position)
{
  bool orientation = track_->orientation();
  if (orientation == SliderTrack::horizontal)
    rel_value_ = (new_position.x() - track_->offset_x() - SliderTrack::stop_offset) /
      (track_->length() - 2.0 * SliderTrack::stop_offset);
  else
  {
    rel_value_ = (track_->offset_y() + track_->length() - SliderTrack::stop_offset - new_position.y()) /
      (track_->length() - 2.0 * SliderTrack::stop_offset);
  }
  rel_value_ = std::clamp(rel_value_, 0.0, 1.0);
  redraw();
}

void Slider::redraw()
{
  Layer* layer = track_->layer();
  handle_ = std::make_shared<draw::SliderHandle>(this);
  highlighted_track_  = std::make_shared<LayerRegion>(DrawHighlightedTrack{this});
  set_value_text();

  layer->draw(value_text_);
  layer->draw(highlighted_track_);
  layer->draw(handle_);
}

void Slider::set_value(double value, double min_value, double max_value)
{
  min_value_ = min_value;
  max_value_ = max_value;
  rel_value_ = std::clamp((value - min_value) / (max_value - min_value), 0.0, 1.0);
  Layer* layer = handle_->layer();
  Window* window = layer->window();
  double pixel_x = track_->offset_x();
  double pixel_y = track_->offset_y();
  if (track_->orientation() == SliderTrack::horizontal)
    pixel_x += horizontal_position_pixels();
  else
    pixel_y += vertical_position_pixels();
  window->update_grabbed(index_, pixel_x, pixel_y);
  redraw();
}

} // namespace cairowindow::draw
