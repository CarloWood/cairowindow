#pragma once

#include "Sample.h"
#include "cairowindow/Point.h"
#include "cairowindow/Text.h"
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif
using utils::has_print_on::operator<<;

// A Sample including plot objects: a point and a label.
class PlotSample
{
 private:
  gradient_descent::Sample const* master_{};

  mutable cairowindow::plot::Point P_;
  mutable cairowindow::plot::Text P_label_;

 public:
  // Required for History.
  PlotSample() = default;

  // Constructor used by LocalExtreme.
  PlotSample(gradient_descent::Sample const* master, cairowindow::plot::Point const& point, cairowindow::plot::Text const& label) :
    master_(master), P_(point), P_label_(label) { }

  // Required for History.
  void initialize(gradient_descent::Sample const* master, cairowindow::plot::Point const& point, cairowindow::plot::Text const& label)
  {
    master_ = master;
    P_ = point;
    P_label_ = label;
  }

  cairowindow::Point const& P() const { return P_; }
  cairowindow::Text const& label() const { return P_label_; }

  gradient_descent::Sample const* sample() const { return master_; }

  double w() const { return master_->w(); }
  double Lw() const { return master_->Lw(); }
  double dLdw() const { return master_->dLdw(); }

#ifdef CWDEBUG
  std::string debug_label() const { return P_label_.text(); }

  void print_on(std::ostream& os) const
  {
    os << debug_label() << " (at w = " << master_->w() << ")";
  }
#endif
};
