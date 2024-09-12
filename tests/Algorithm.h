#pragma once

#include "gradient_descent2/Algorithm.h"
#include "EnableDrawing.h"
#include "cwds/Restart.h"
#include <memory>

namespace enable_drawing {

class Algorithm : public gradient_descent::Algorithm
{
 private:
  std::unique_ptr<EnableDrawing> enable_drawing_;

 public:
  using gradient_descent::Algorithm::Algorithm;

  void enable_drawing(Function const& L, double w_min, double w_max)
  {
    enable_drawing_ = std::make_unique<EnableDrawing>(this, L, w_min, w_max);
  }

  bool operator()(double& w, double Lw, double dLdw)
  {
    bool result = gradient_descent::Algorithm::operator()(w, Lw, dLdw);
    if (enable_drawing_ && !debug::Restart<0>::s_restarting)
      enable_drawing_->wait();
    return result;
  }

  EnableDrawing& enable_drawing() { return *enable_drawing_; }
  EnableDrawing const& enable_drawing() const { return *enable_drawing_; }
};

} // namespace enable_drawing