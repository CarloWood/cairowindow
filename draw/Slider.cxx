#include "sys.h"
#include "Slider.h"
#include "cairowindow/Layer.h"
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

void Slider::draw_regions_on(Layer* layer)
{
  DoutEntering(dc::cairowindow, "Slider::draw_regions_on(" << layer << ")");

  highlighted_track_ = std::make_shared<LayerRegion>(DrawHighlightedTrack{this});

  layer->draw(track_);
  layer->draw(highlighted_track_);
  layer->draw(handle_);
}

void Slider::moved(plot::Plot* plot, cairowindow::Point const& new_position)
{
  bool orientation = track_->orientation();
  if (orientation == SliderTrack::horizontal)
    value_ = (new_position.x() - track_->offset_x() - SliderTrack::stop_offset) /
      (track_->length() - 2.0 * SliderTrack::stop_offset);
  else
  {
    value_ = (track_->offset_y() + track_->length() - SliderTrack::stop_offset - new_position.y()) /
      (track_->length() - 2.0 * SliderTrack::stop_offset);
  }
  value_ = std::clamp(value_, 0.0, 1.0);

  Layer* layer = track_->layer();
  handle_ = std::make_shared<draw::SliderHandle>(this);
  highlighted_track_  = std::make_shared<LayerRegion>(DrawHighlightedTrack{this});
  layer->draw(highlighted_track_);
  layer->draw(handle_);
}

} // namespace cairowindow::draw
