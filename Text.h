#pragma once

#include "Pixel.h"
#include "utils/Badge.h"
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

template<CS> class CoordinateSystem;
template<CS> class CoordinateMapper;

namespace plot {
class Plot;

class Text : public cairowindow::Text
{
 public:
  explicit Text(cairowindow::Text const& text) : cairowindow::Text(text) { }
  using cairowindow::Text::Text;

 public:
  template<CS cs, typename... Args>
  void create_draw_object(utils::Badge<Plot, cairowindow::CoordinateSystem<cs>, cairowindow::CoordinateMapper<cs>>, Args&&... args) const
  {
    draw_object_ = std::make_shared<cairowindow::draw::Text>(std::forward<Args>(args)...);
  }

  // Accessor for the draw object; used by Plot and CoordinateSystem.
  std::shared_ptr<draw::Text> const& draw_object() const
  {
    return draw_object_;
  }

  void reset()
  {
    draw_object_.reset();
  }

 private:
  mutable std::shared_ptr<draw::Text> draw_object_;
};

} // namespace plot
} // namespace cairowindow
