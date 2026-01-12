#pragma once

#include "Line.h"
#include "Direction.h"
#include <memory>

namespace cairowindow {

class Arc
{
 public:
  using Line = math::Line<2, double>;

 private:
  Point center_;
  double start_angle_ = 0.0;
  double end_angle_ = 0.0;
  double radius_ = 0.0;

 public:
  Arc(Point const& center, double start_angle, double end_angle, double radius = 0.0) :
    center_(center), start_angle_(start_angle), end_angle_(end_angle), radius_(radius) { }

  Arc(Point const& center, Direction const& start, Direction const& end, double radius = 0.0) :
    center_(center), start_angle_(start.as_angle()), end_angle_(end.as_angle()), radius_(radius) { }

  Arc(Line const& line1, Line const& line2, double radius = 0.0) :
    center_(line1.intersection_with(line2)), start_angle_(line1.direction().as_angle()),
    end_angle_(line2.direction().as_angle()), radius_(radius) { }

  bool is_defined() const { return start_angle_ != end_angle_; }
  bool has_radius() const { return radius_ != 0.0; }

  Point const& center() const { return center_; }
  double start_angle() const { return start_angle_; }
  double end_angle() const { return end_angle_; }
  double radius() const { return radius_; }

  Direction bisector_direction() const
  {
    double angle_diff = end_angle_ - start_angle_;
    if (angle_diff < 0.0)
      angle_diff += 2 * M_PI;
    return Direction{start_angle_ + 0.5 * angle_diff};
  }
};

namespace draw {
class Arc;
} // namespace draw

namespace plot {
class Plot;

class Arc : public cairowindow::Arc
{
 public:
  explicit Arc(cairowindow::Arc const& arc) : cairowindow::Arc(arc) { }
  using cairowindow::Arc::Arc;

 public:
  friend class Plot;
  mutable std::shared_ptr<draw::Arc> draw_object_;
};

} // namespace plot
} // namespace cairowindow
