#pragma once

namespace cairowindow::plot {

class Range
{
 private:
  double min_;
  double max_;

 public:
  Range() : min_{0.0}, max_{1.0} { }
  Range(double min, double max) : min_(min), max_(max) { }

  double min() const { return min_; }
  double max() const { return max_; }
};

} // namespace cairowindow::plot
