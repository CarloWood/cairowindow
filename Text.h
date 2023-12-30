#pragma once

namespace cairowindow {

class Text
{
 private:
  Point position_;
  std::string text_;

 public:
  Text(Point const& position, std::string const& text) : position_(position), text_(text) { }

  Point const& position() const { return position_; }
  std::string const& text() const { return text_; }
};

} // namespace cairowindow
