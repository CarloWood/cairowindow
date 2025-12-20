#pragma once

#include "Shape.h"
#include "Text.h"
#include "MultiRegion.h"
#include "cairowindow/Point.h"
#include "cairowindow/Color.h"
#include "cairowindow/Draggable.h"
#include "cairowindow/debug_channel.h"

namespace cairowindow::draw {

class Slider;
struct DrawHighlightedTrack;

class SliderTrack : public LayerRegion
{
 public:
  static constexpr bool horizontal = false;
  static constexpr bool vertical = true;
  static constexpr int bg_width = 32;
  static constexpr int track_width = 5;
  static constexpr int stop_offset = 15;

 private:
  bool orientation_;
  double offset_x_;             // The coordinates of the middle point on one side of the slider.
  double offset_y_;
  double width_;                // The width of the knob.
  double length_;
  //
  //             .-----.-.------------------------.
  //  offset --> + ====| |======================= | > track_width
  //             '-----'-'------------------------'
  //              <------------length------------>

 public:
  SliderTrack(double offset_x, double offset_y, double width, double height) :
    orientation_((width < height) ? vertical : horizontal),
    offset_x_(offset_x), offset_y_(offset_y), width_(std::min(width, height)), length_(std::max(width, height))
  {
  }

  void draw_highlighted_track(cairo_t* cr, Slider const* parent) const;

  bool orientation() const { return orientation_; }
  double offset_x() const { return offset_x_; }
  double offset_y() const { return offset_y_; }
  double width() const { return width_; }
  double length() const { return length_; }

 private:
  StrokeExtents do_draw(cairo_t* cr) override
  {
    DoutEntering(dc::cairowindow, "draw::Slider::do_draw(cr) [" << this << "]");
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    constexpr Color bg_color(0.498, 0.498, 0.498);
    constexpr Color track_color1(0.2627, 0.2627, 0.2627);
    constexpr Color border_color(0.0588, 0.0588, 0.0588);

    cairo_set_source_rgb(cr, bg_color.red(), bg_color.green(), bg_color.blue());
    if (orientation_ == vertical)
      cairo_rectangle(cr, offset_x_ - 0.5 * bg_width, offset_y_, bg_width, length_);
    else
      cairo_rectangle(cr, offset_x_, offset_y_ - 0.5 * bg_width, length_, bg_width);

    double x1, y1, x2, y2;
    cairo_fill_extents(cr, &x1, &y1, &x2, &y2);
    cairo_fill(cr);

    cairo_set_source_rgb(cr, track_color1.red(), track_color1.green(), track_color1.blue());
    if (orientation_ == vertical)
      cairo_rectangle(cr, offset_x_ - 0.5 * track_width, offset_y_ + stop_offset, track_width, length_ - 2 * stop_offset);
    else
      cairo_rectangle(cr, offset_x_ + stop_offset, offset_y_ - 0.5 * track_width, length_ - 2 * stop_offset, track_width);

    cairo_fill(cr);

    return {x1, y1, x2, y2};
  }
};

class SliderHandle : public Shape
{
 public:
  SliderHandle(Slider const* parent);
};

class Slider : public MultiRegion, public plot::Draggable
{
 public:
  static constexpr Color highlighted_color{0.098, 0.596, 0.596};

 private:
  double rel_value_;                                    // Slider value in the range of [0, 1].
  double min_value_;                                    // The value corresponding to a rel_value of 0.
  double max_value_;                                    // Idem 1.
  std::shared_ptr<SliderTrack> track_;
  std::shared_ptr<Text> value_text_;
  std::shared_ptr<LayerRegion> highlighted_track_;
  std::shared_ptr<SliderHandle> handle_;

 public:
  Slider(double x, double y, double width, double height, double start_value, double min_value, double max_value) :
    MultiRegion(highlighted_color, 0.0),
    rel_value_(std::clamp((start_value - min_value) / (max_value - min_value), 0.0, 1.0)), min_value_(min_value), max_value_(max_value),
    track_(std::make_shared<SliderTrack>(x, y, width, height)), handle_(std::make_shared<SliderHandle>(this)) { }

  SliderTrack const& track() const { return *track_; }
  SliderHandle const& handle() const { return *handle_; }
  double rel_value() const { return rel_value_; }
  double horizontal_position_pixels() const
  {
    return SliderTrack::stop_offset + rel_value_ * (track_->length() - 2.0 * SliderTrack::stop_offset);
  }
  double vertical_position_pixels() const
  {
    return track_->length() - SliderTrack::stop_offset - rel_value_ * (track_->length() - 2.0 * SliderTrack::stop_offset);
  }

  // Adjust values.
  void set_value(double value, double min_value, double max_value);

 private:
  friend struct DrawHighlightedTrack;
  void do_set_fill_color(cairo_t* cr) const { return set_fill_color(cr); }
  StrokeExtents do_fill(cairo_t* cr) const { return fill(cr); }

 private:
  void set_value_text();
  void redraw();

  // Implementation of MultiRegion.
  void draw_regions_on(Layer* layer) override;

  // Implementation of plot::Draggable.
  cairowindow::Rectangle const& geometry() const override { return handle_->geometry(); }
  void moved(plot::Plot* plot, cairowindow::Point const& new_position) override;
  void set_position(cairowindow::Point const& new_position) override
  {
    // Moving a slider with a Point value not implemented.
    ASSERT(false);
  }
  bool convert() const override { return false; }       // We want mouse coordinates.

#ifdef CWDEBUG
  void print_on(std::ostream& os) const override
  {
    os << "draw::Slider{" << (100.0 * rel_value_) << "%}";
  }
#endif
};

} // namespace cairowindow::draw
