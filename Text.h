#pragma once

#include "Pixel.h"
#include <string>

namespace cairowindow {

class Text
{
 private:
  Pixel position_;
  std::string text_;

 public:
  Text(Pixel position, std::string const& text) : position_(position), text_(text) { }

  Pixel position() const { return position_; }
  std::string const& text() const { return text_; }
};

} // namespace cairowindow
