#pragma once

#include "Point.h"
#include "BezierFitter.h"
#include <vector>

namespace cairowindow {

class Curve
{
 private:
  std::vector<Point> points_;

 public:
  Curve() = default;
  Curve(std::vector<Point> const& points) : points_(points) { }
  Curve(std::vector<Point>&& points) : points_(std::move(points)) { }

  std::vector<Point> const& points() const { return points_; }
  std::vector<Point>& points() { return points_; }
};

class FunctionCurve
{
 private:
  BezierFitter fitter_;
};

} // namespace cairowindow
