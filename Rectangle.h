#pragma once

#include "debug.h"

namespace cairowindow {

class Rectangle
{
 private:
  double offset_x_;
  double offset_y_;
  double width_;
  double height_;

 public:
  Rectangle() : offset_x_(0.0), offset_y_(0.0), width_(0.0), height_(0.0) { }
  Rectangle(double offset_x, double offset_y, double width, double height) :
    offset_x_(offset_x), offset_y_(offset_y), width_(width), height_(height) { }

  double offset_x() const { return offset_x_; }
  double offset_y() const { return offset_y_; }
  double width() const { return width_; }
  double height() const { return height_; }

  bool is_defined() const { return width_ > 0.0 && height_ > 0.0; }

  double area() const { ASSERT(is_defined()); return width_ * height_; }

  bool overlaps(Rectangle const& rectangle) const
  {
    double a_left = offset_x_;
    double a_right = offset_x_ + width_;
    double a_top = offset_y_ + height_;
    double a_bottom = offset_y_;

    double b_left = rectangle.offset_x_;
    double b_right = rectangle.offset_x_ + rectangle.width_;
    double b_top = rectangle.offset_y_ + rectangle.height_;
    double b_bottom = rectangle.offset_y_;

    if (a_right < b_left || b_right < a_left) return false;
    if (a_top < b_bottom || b_top < a_bottom) return false;
    return true;
  }
};

} // namespace cairowindow
