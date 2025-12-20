#pragma once

#include "Point.h"
#include <memory>

namespace cairowindow {

class Circle
{
 private:
  Point center_;
  double radius_;

 public:
  Circle() = default;
  Circle(Point const& center, double radius) : center_(center), radius_(radius) { }

  Point const& center() const { return center_; }
  double radius() const { return radius_; }
};

namespace draw {
class Circle;
} // namespace draw

namespace plot {
class Plot;

//--------------------------------------------------------------------------
// Circle

class Circle : public cairowindow::Circle
{
 public:
  explicit Circle(cairowindow::Circle const& circle) : cairowindow::Circle(circle) { }
  using cairowindow::Circle::Circle;

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Circle> draw_object_;
};

} // namespace plot
} // namespace cairowindow

