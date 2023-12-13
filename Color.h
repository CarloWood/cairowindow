#pragma once

namespace cairowindow {

class Color
{
 private:
  float red_;
  float green_;
  float blue_;
  float alpha_;

 public:
  Color(float red, float green, float blue, float alpha = 1.0f) : red_(red), green_(green), blue_(blue), alpha_(alpha) { }

  float red() const { return red_; }
  float green() const { return green_; }
  float blue() const { return blue_; }
  float alpha() const { return alpha_; }

  bool is_opaque() const { return alpha_ == 1.0f; }
};

} // namespace cairowindow
