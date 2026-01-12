#pragma once

#include "Rectangle.h"
#include "Color.h"

namespace cairowindow {

class LayerArgs
{
 protected:
  Geometry geometry_;

 public:
  LayerArgs() { }
  LayerArgs(Geometry geometry) : geometry_(geometry) { }

  bool has_geometry() const { return geometry_.is_defined(); }
  Geometry const& geometry() const { return geometry_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{geometry_:" << geometry_ << '}';
  }
#endif
};

class BackgroundLayerArgs : public LayerArgs
{
 private:
  Color background_color_;

 public:
  BackgroundLayerArgs(Color background_color) : background_color_(background_color) { }
  BackgroundLayerArgs(Geometry geometry, Color background_color) : LayerArgs(geometry), background_color_(background_color) { }

  Color const& background_color() const { return background_color_; }

#ifdef CWDEBUG
  void print_on(std::ostream& os) const
  {
    os << "{background_color_:" << background_color_ << ", geometry_:" << geometry_ << '}';
  }
#endif
};

} // namespace cairowindow
