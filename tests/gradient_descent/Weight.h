#pragma once

namespace gradient_descent {

// This class represents a point for which we consider
// to calculate the next L(w) and L'(w) values.
//
class Weight
{
 protected:
  double w_;                    // w

 public:
  // Required by Sample() and find_extreme.
  Weight() = default;

  Weight(double w) : w_(w) { }

  Weight& operator-=(double delta_w)
  {
    w_ -= delta_w;
    return *this;
  }

  Weight& operator+=(double delta_w)
  {
    w_ += delta_w;
    return *this;
  }

  operator double() const
  {
    return w_;
  }
};

} // namespace gradient_descent
