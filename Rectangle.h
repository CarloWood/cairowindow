#pragma once

namespace cairowindow {

class Rectangle
{
 private:
  double offset_x_;
  double offset_y_;
  double width_;
  double height_;

 public:
  Rectangle(double offset_x, double offset_y, double width, double height) :
    offset_x_(offset_x), offset_y_(offset_y), width_(width), height_(height) { }

  double offset_x() const { return offset_x_; }
  double offset_y() const { return offset_y_; }
  double width() const { return width_; }
  double height() const { return height_; }
};

} // namespace cairowindow
