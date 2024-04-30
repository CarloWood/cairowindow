#pragma once

#include "Pixel.h"
#include <string>
#include <memory>

namespace cairowindow {

class Text
{
 private:
  Pixel position_;
  std::string text_;

 public:
  Text() = default;
  Text(Pixel position, std::string const& text) : position_(position), text_(text) { }

  Pixel position() const { return position_; }
  std::string const& text() const { return text_; }
};

namespace draw {
class Text;
} // namespace draw

namespace plot {
class Plot;

class Text : public cairowindow::Text
{
 public:
  using cairowindow::Text::Text;
  Text(cairowindow::Text const& text) : cairowindow::Text(text) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Text> draw_object_;
};

} // namespace plot
} // namespace cairowindow
