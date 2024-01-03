#pragma once

namespace cairowindow {

class Direction;

class Point
{
 private:
  double x_;
  double y_;

 public:
  Point(double x, double y) : x_(x), y_(y) { }

  double x() const { return x_; }
  double y() const { return y_; }

  Point operator+(Direction const& direction);
};

} // namespace cairowindow
