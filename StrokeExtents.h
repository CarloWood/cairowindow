#pragma once

#include "Rectangle.h"
#include "utils/has_print_on.h"
#include <cairo/cairo.h>
#include <cmath>

namespace cairowindow {
using utils::has_print_on::operator<<;

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
  StrokeExtents(double x1, double y1, double x2, double y2) :
    half_width_(0.5 * (x2 - x1)), half_height_(0.5 * (y2 - y1)), center_x_(0.5 * (x2 + x1)), center_y_(0.5 * (y2 + y1)) { }
  StrokeExtents(Rectangle const& rectangle) : half_width_(0.5 * rectangle.width()), half_height_(0.5 * rectangle.height()),
    center_x_(rectangle.offset_x() + half_width_), center_y_(rectangle.offset_y() + half_height_) { }

  bool is_defined() const { return half_width_ != 0.0 || half_height_ != 0.0; }
  double area() const { return 4.0 * half_width_ * half_height_; }

  void set_path(cairo_t* cr) const
  {
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

} // namespace cairowindow
