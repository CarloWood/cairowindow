#pragma once

#include "Rectangle.h"
#include "Color.h"

namespace cairowindow {

class LayerArgs
{
 private:
  Rectangle rectangle_;

 public:
  LayerArgs() { }
  LayerArgs(Rectangle rectangle) : rectangle_(rectangle) { }

  bool has_rectangle() const { return rectangle_.is_defined(); }
  Rectangle const& rectangle() const { return rectangle_; }
};

class BackgroundLayerArgs : public LayerArgs
{
 private:
  Color background_color_;

 public:
  BackgroundLayerArgs(Color background_color) : background_color_(background_color) { }
  BackgroundLayerArgs(Rectangle rectangle, Color background_color) : LayerArgs(rectangle), background_color_(background_color) { }

  Color const& background_color() const { return background_color_; }
};

} // namespace cairowindow
