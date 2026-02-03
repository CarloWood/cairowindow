#pragma once

#include "IntersectRectangle.h"
#include <cairo/cairo.h>
#include <cmath>
#include "debug.h"
#ifdef CWDEBUG
#include "cairowindow/debugcairo.h"
#include "utils/has_print_on.h"
#endif

namespace cairowindow {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class StrokeExtents
{
 private:
  double half_width_;
  double half_height_;
  double center_x_;
  double center_y_;

 public:
  // Construct an "undefined" StrokeExtents object (with half_width_ and half_height_ zero).
  StrokeExtents() : half_width_(0.0), half_height_(0.0) { }
  StrokeExtents(double x1, double y1, double x2, double y2)
  {
    // Include anti-aliasing pixels in extent.
    x1 = std::floor(x1);
    y1 = std::floor(y1);
    x2 = std::ceil(x2);
    y2 = std::ceil(y2);

    // It is acceptable that a StrokeExtents falls outside of the window;
    // this will be tested upon return. However, it is not ok to have reversed
    // x or y coordinates.
    ASSERT(x2 >= x1 && y2 >= y1);

    half_width_ = 0.5 * (x2 - x1);
    half_height_ = 0.5 * (y2 - y1);
    center_x_ = 0.5 * (x2 + x1);
    center_y_ = 0.5 * (y2 + y1);
  }
  StrokeExtents(Geometry const& geometry) : StrokeExtents(geometry.offset_x(), geometry.offset_y(),
      geometry.offset_x() + geometry.width(), geometry.offset_y() + geometry.height()) { }

  bool clip(Geometry const& geometry)
  {
    IntersectRectangle<csid::pixels> intersection(geometry, *this);
    double width = intersection.x2() - intersection.x1();
    double height = intersection.y2() - intersection.y1();
    if (width <= 0.0 || height <= 0.0)
      return false;
    half_width_ = 0.5 * width;
    half_height_ = 0.5 * height;
    center_x_ = 0.5 * (intersection.x1() + intersection.x2());
    center_y_ = 0.5 * (intersection.y1() + intersection.y2());
    return true;
  }

  // Accessors.
  double center_x() const { return center_x_; }
  double center_y() const { return center_y_; }
  double width() const { return 2.0 * half_width_; }
  double height() const { return 2.0 * half_height_; }
  double x1() const { return center_x_ - half_width_; }
  double y1() const { return center_y_ - half_height_; }

  bool is_defined() const { return half_width_ != 0.0 || half_height_ != 0.0; }
  double area() const { return 4.0 * half_width_ * half_height_; }

  void set_path(cairo_t* cr) const
  {
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    cairo_rectangle(cr, center_x_ - half_width_, center_y_ - half_height_, 2.0 * half_width_, 2.0 * half_height_);
  }

  void unpack(int& x_out, int& y_out, int& width_out, int& height_out) const
  {
    x_out = static_cast<int>(std::floor(center_x_ - half_width_));
    y_out = static_cast<int>(std::floor(center_y_ - half_height_));
    int xm = static_cast<int>(std::ceil(center_x_ + half_width_));
    int ym = static_cast<int>(std::ceil(center_y_ + half_height_));
    width_out = xm - x_out;
    height_out = ym - y_out;
  }

  void print_on(std::ostream& os) const
  {
    os << '{' << "half_width_:" << half_width_ << ", half_height_:" << half_height_ <<
      ", center_x_:" << center_x_ << ", center_y_:" << center_y_ << '}';
  }

  bool overlaps(StrokeExtents const& stroke_extents) const
  {
    return std::abs(center_x_ - stroke_extents.center_x_) < half_width_ + stroke_extents.half_width_ &&
           std::abs(center_y_ - stroke_extents.center_y_) < half_height_ + stroke_extents.half_height_;
  }
};

// Specialization for csid::pixels.
template<CS cs>
requires (cs == csid::pixels || cs == csid::plot)
IntersectRectangle<cs>::IntersectRectangle(StrokeExtents const& stroke_extents) requires (cs == csid::pixels) :
  x1_(stroke_extents.x1()), y1_(stroke_extents.y1()),
  x2_(x1_ + stroke_extents.width()), y2_(y1_ + stroke_extents.height()) { }

} // namespace cairowindow
