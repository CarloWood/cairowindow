#pragma once

namespace cairowindow {

class Color
{
 private:
  double red_;
  double green_;
  double blue_;
  double alpha_;

 public:
  constexpr Color(double red, double green, double blue, double alpha = 1.0f) : red_(red), green_(green), blue_(blue), alpha_(alpha) { }

  double red() const { return red_; }
  double green() const { return green_; }
  double blue() const { return blue_; }
  double alpha() const { return alpha_; }

  bool is_opaque() const { return alpha_ == 1.0f; }
};

namespace color {

constexpr Color red{255, 0, 0};
constexpr Color green{0, 255, 0};
constexpr Color blue{0, 0, 255};
constexpr Color yellow{255, 255, 0};
constexpr Color magenta{255, 0, 255};
constexpr Color cyan{0, 255, 255};
constexpr Color black{0, 0, 0};
constexpr Color white{255, 255, 255};
constexpr Color gray{128, 128, 128};
constexpr Color orange{255, 165, 0};
constexpr Color brown{165, 42, 42};
constexpr Color pink{255, 192, 203};
constexpr Color violet{238, 130, 238};
constexpr Color indigo{75, 0, 130};
constexpr Color maroon{128, 0, 0};
constexpr Color turquoise{64, 224, 208};
constexpr Color tan{210, 180, 140};
constexpr Color gold{255, 215, 0};
constexpr Color beige{245, 245, 220};
constexpr Color lavender{230, 230, 250};
constexpr Color peach{255, 218, 185};
constexpr Color lime{0, 255, 0};
constexpr Color olive{128, 128, 0};
constexpr Color teal{0, 128, 128};
constexpr Color navy{0, 0, 128};
constexpr Color coral{255, 127, 80};
constexpr Color salmon{250, 128, 114};

} // namespace color

} // namespace cairowindow
