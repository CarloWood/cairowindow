#pragma once

#include "math/cs/Line.h"
#include <cmath>

namespace cairowindow::cs {

template<CS cs>
class Arc
{
 public:
  using Line = math::cs::Line<cs>;

 private:
  math::cs::Point<cs> center_;
  double start_angle_ = 0.0;
  double end_angle_ = 0.0;
  double radius_ = 0.0;

 public:
  // Construct an undefined arc.
  Arc() = default;

  Arc(math::cs::Point<cs> const& center, double start_angle, double end_angle, double radius = 0.0) :
    center_(center), start_angle_(start_angle), end_angle_(end_angle), radius_(radius) { }

  Arc(math::cs::Point<cs> const& center, math::cs::Direction<cs> const& start, math::cs::Direction<cs> const& end, double radius = 0.0) :
    center_(center), start_angle_(start.as_angle()), end_angle_(end.as_angle()), radius_(radius) { }

  Arc(math::cs::Line<cs> const& line1, math::cs::Line<cs> const& line2, double radius = 0.0) :
    center_(line1.intersection_with(line2)), start_angle_(line1.direction().as_angle()),
    end_angle_(line2.direction().as_angle()), radius_(radius) { }

  bool is_defined() const { return start_angle_ != end_angle_; }
  bool has_radius() const { return radius_ != 0.0; }

  math::cs::Point<cs> const& center() const { return center_; }
  double start_angle() const { return start_angle_; }
  double end_angle() const { return end_angle_; }
  double radius() const { return radius_; }

  math::cs::Direction<cs> bisector_direction() const
  {
    double angle_diff = end_angle_ - start_angle_;
    if (angle_diff < 0.0)
      angle_diff += 2 * M_PI;
    return math::cs::Direction<cs>{start_angle_ + 0.5 * angle_diff};
  }
};

} // namespace cairowindow::cs
