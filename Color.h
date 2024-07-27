#pragma once

#ifdef CWDEBUG
#include "utils/has_print_on.h"
#include <iostream>
#endif
#include "debug.h"

namespace cairowindow {
#ifdef CWDEBUG
// This class defines a print_on method.
using utils::has_print_on::operator<<;
#endif

class Color
{
 private:
  double red_;
  double green_;
  double blue_;
  double alpha_;

 public:
  // Construct an "undefined" color.
  constexpr Color() : red_(2.0), green_(2.0), blue_(2.0), alpha_(2.0) { }
  constexpr Color(double red, double green, double blue, double alpha = 1.0) : red_(red), green_(green), blue_(blue), alpha_(alpha) { }

  double red() const { return red_; }
  double green() const { return green_; }
  double blue() const { return blue_; }
  double alpha() const { return alpha_; }

  bool is_defined() const { return alpha_ != 2.0; }
  bool is_transparent() const { return alpha_ == 0.0; }
  bool is_opaque() const { return alpha_ == 1.0; }
  void set_opaque() { alpha_ = 1.0; }

  // Get the (next) color from a color pool.
  static Color get_color(int color_index);
  static Color next_color();

  // The Style macros need an operator!= to compare with the "undefined" magic.
  bool operator!=(Color const& undefined_color) const
  {
    // Only use this to compare with a default constructed color.
    ASSERT(!undefined_color.is_defined());
    return is_defined();
  }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{red_:" << red_ << ", green_:" << green_ << ", blue_:" << blue_ << ", alpha_:" << alpha_ << '}';
  }
#endif
};

namespace color {

constexpr Color undefined;

// Colors known to Qt, copied from https://www.w3.org/TR/SVG11/types.html#ColorKeywords .
constexpr Color aliceblue{0.941, 0.973, 1.0};
constexpr Color antiquewhite{0.98, 0.922, 0.843};
constexpr Color aqua{0.0, 1.0, 1.0};
constexpr Color aquamarine{0.498, 1.0, 0.831};
constexpr Color azure{0.941, 1.0, 1.0};
constexpr Color beige{0.961, 0.961, 0.863};
constexpr Color bisque{1.0, 0.894, 0.769};
constexpr Color black{0.0, 0.0, 0.0};
constexpr Color blanchedalmond{1.0, 0.922, 0.804};
constexpr Color blue{0.0, 0.0, 1.0};
constexpr Color blueviolet{0.541, 0.169, 0.886};
constexpr Color brown{0.647, 0.165, 0.165};
constexpr Color burlywood{0.871, 0.722, 0.529};
constexpr Color cadetblue{0.373, 0.62, 0.627};
constexpr Color chartreuse{0.498, 1.0, 0.0};
constexpr Color chocolate{0.824, 0.412, 0.118};
constexpr Color coral{1.0, 0.498, 0.314};
constexpr Color cornflowerblue{0.392, 0.584, 0.929};
constexpr Color cornsilk{1.0, 0.973, 0.863};
constexpr Color crimson{0.863, 0.078, 0.235};
constexpr Color cyan{0.0, 1.0, 1.0};
constexpr Color darkblue{0.0, 0.0, 0.545};
constexpr Color darkcyan{0.0, 0.545, 0.545};
constexpr Color darkgoldenrod{0.722, 0.525, 0.043};
constexpr Color darkgray{0.663, 0.663, 0.663};
constexpr Color darkgreen{0.0, 0.392, 0.0};
constexpr Color darkgrey{0.663, 0.663, 0.663};
constexpr Color darkkhaki{0.741, 0.718, 0.42};
constexpr Color darkmagenta{0.545, 0.0, 0.545};
constexpr Color darkolivegreen{0.333, 0.42, 0.184};
constexpr Color darkorange{1.0, 0.549, 0.0};
constexpr Color darkorchid{0.6, 0.196, 0.8};
constexpr Color darkred{0.545, 0.0, 0.0};
constexpr Color darksalmon{0.914, 0.588, 0.478};
constexpr Color darkseagreen{0.561, 0.737, 0.561};
constexpr Color darkslateblue{0.282, 0.239, 0.545};
constexpr Color darkslategray{0.184, 0.31, 0.31};
constexpr Color darkslategrey{0.184, 0.31, 0.31};
constexpr Color darkturquoise{0.0, 0.808, 0.82};
constexpr Color darkviolet{0.58, 0.0, 0.827};
constexpr Color deeppink{1.0, 0.078, 0.576};
constexpr Color deepskyblue{0.0, 0.749, 1.0};
constexpr Color dimgray{0.412, 0.412, 0.412};
constexpr Color dimgrey{0.412, 0.412, 0.412};
constexpr Color dodgerblue{0.118, 0.565, 1.0};
constexpr Color firebrick{0.698, 0.133, 0.133};
constexpr Color floralwhite{1.0, 0.98, 0.941};
constexpr Color forestgreen{0.133, 0.545, 0.133};
constexpr Color fuchsia{1.0, 0.0, 1.0};
constexpr Color gainsboro{0.863, 0.863, 0.863};
constexpr Color ghostwhite{0.973, 0.973, 1.0};
constexpr Color gold{1.0, 0.843, 0.0};
constexpr Color goldenrod{0.855, 0.647, 0.125};
constexpr Color gray{0.502, 0.502, 0.502};
constexpr Color grey{0.502, 0.502, 0.502};
constexpr Color green{0.0, 0.502, 0.0};
constexpr Color greenyellow{0.678, 1.0, 0.184};
constexpr Color honeydew{0.941, 1.0, 0.941};
constexpr Color hotpink{1.0, 0.412, 0.706};
constexpr Color indianred{0.804, 0.361, 0.361};
constexpr Color indigo{0.294, 0.0, 0.51};
constexpr Color ivory{1.0, 1.0, 0.941};
constexpr Color khaki{0.941, 0.902, 0.549};
constexpr Color lavender{0.902, 0.902, 0.98};
constexpr Color lavenderblush{1.0, 0.941, 0.961};
constexpr Color lawngreen{0.486, 0.988, 0.0};
constexpr Color lemonchiffon{1.0, 0.98, 0.804};
constexpr Color lightblue{0.678, 0.847, 0.902};
constexpr Color lightcoral{0.941, 0.502, 0.502};
constexpr Color lightcyan{0.878, 1.0, 1.0};
constexpr Color lightgoldenrodyellow{0.98, 0.98, 0.824};
constexpr Color lightgray{0.827, 0.827, 0.827};
constexpr Color lightgreen{0.565, 0.933, 0.565};
constexpr Color lightgrey{0.827, 0.827, 0.827};
constexpr Color lightpink{1.0, 0.714, 0.757};
constexpr Color lightred{1.0, 0.5, 0.5};                // Not in the database.
constexpr Color lightsalmon{1.0, 0.627, 0.478};
constexpr Color lightseagreen{0.125, 0.698, 0.667};
constexpr Color lightskyblue{0.529, 0.808, 0.98};
constexpr Color lightslategray{0.467, 0.533, 0.6};
constexpr Color lightslategrey{0.467, 0.533, 0.6};
constexpr Color lightsteelblue{0.69, 0.769, 0.871};
constexpr Color lightyellow{1.0, 1.0, 0.878};
constexpr Color lime{0.0, 1.0, 0.0};
constexpr Color limegreen{0.196, 0.804, 0.196};
constexpr Color linen{0.98, 0.941, 0.902};
constexpr Color magenta{1.0, 0.0, 1.0};
constexpr Color maroon{0.502, 0.0, 0.0};
constexpr Color mediumaquamarine{0.4, 0.804, 0.667};
constexpr Color mediumblue{0.0, 0.0, 0.804};
constexpr Color mediumorchid{0.729, 0.333, 0.827};
constexpr Color mediumpurple{0.576, 0.439, 0.859};
constexpr Color mediumseagreen{0.235, 0.702, 0.443};
constexpr Color mediumslateblue{0.482, 0.408, 0.933};
constexpr Color mediumspringgreen{0.0, 0.98, 0.604};
constexpr Color mediumturquoise{0.282, 0.82, 0.8};
constexpr Color mediumvioletred{0.78, 0.082, 0.522};
constexpr Color midnightblue{0.098, 0.098, 0.439};
constexpr Color mintcream{0.961, 1.0, 0.98};
constexpr Color mistyrose{1.0, 0.894, 0.882};
constexpr Color moccasin{1.0, 0.894, 0.71};
constexpr Color navajowhite{1.0, 0.871, 0.678};
constexpr Color navy{0.0, 0.0, 0.502};
constexpr Color oldlace{0.992, 0.961, 0.902};
constexpr Color olive{0.502, 0.502, 0.0};
constexpr Color olivedrab{0.42, 0.557, 0.137};
constexpr Color orange{1.0, 0.647, 0.0};
constexpr Color orangered{1.0, 0.271, 0.0};
constexpr Color orchid{0.855, 0.439, 0.839};
constexpr Color palegoldenrod{0.933, 0.91, 0.667};
constexpr Color palegreen{0.596, 0.984, 0.596};
constexpr Color paleturquoise{0.686, 0.933, 0.933};
constexpr Color palevioletred{0.859, 0.439, 0.576};
constexpr Color papayawhip{1.0, 0.937, 0.835};
constexpr Color peachpuff{1.0, 0.855, 0.725};
constexpr Color peru{0.804, 0.522, 0.247};
constexpr Color pink{1.0, 0.753, 0.796};
constexpr Color plum{0.867, 0.627, 0.867};
constexpr Color powderblue{0.69, 0.878, 0.902};
constexpr Color purple{0.502, 0.0, 0.502};
constexpr Color red{1.0, 0.0, 0.0};
constexpr Color rosybrown{0.737, 0.561, 0.561};
constexpr Color royalblue{0.255, 0.412, 0.882};
constexpr Color saddlebrown{0.545, 0.271, 0.075};
constexpr Color salmon{0.98, 0.502, 0.447};
constexpr Color sandybrown{0.957, 0.643, 0.376};
constexpr Color seagreen{0.18, 0.545, 0.341};
constexpr Color seashell{1.0, 0.961, 0.933};
constexpr Color sienna{0.627, 0.322, 0.176};
constexpr Color silver{0.753, 0.753, 0.753};
constexpr Color skyblue{0.529, 0.808, 0.922};
constexpr Color slateblue{0.416, 0.353, 0.804};
constexpr Color slategray{0.439, 0.502, 0.565};
constexpr Color slategrey{0.439, 0.502, 0.565};
constexpr Color snow{1.0, 0.98, 0.98};
constexpr Color springgreen{0.0, 1.0, 0.498};
constexpr Color steelblue{0.275, 0.51, 0.706};
constexpr Color tan{0.824, 0.706, 0.549};
constexpr Color teal{0.0, 0.502, 0.502};
constexpr Color thistle{0.847, 0.749, 0.847};
constexpr Color tomato{1.0, 0.388, 0.278};
constexpr Color turquoise{0.251, 0.878, 0.816};
constexpr Color violet{0.933, 0.51, 0.933};
constexpr Color wheat{0.961, 0.871, 0.702};
constexpr Color white{1.0, 1.0, 1.0};
constexpr Color whitesmoke{0.961, 0.961, 0.961};
constexpr Color yellow{1.0, 1.0, 0.0};
constexpr Color yellowgreen{0.604, 0.804, 0.196};

constexpr Color transparent{0.0, 0.0, 0.0, 0.0};

} // namespace color

} // namespace cairowindow
