#pragma once

#include "MultiRegion.h"
#include "Line.h"

namespace cairowindow::draw {

struct ChessDiagramStyle
{
  Color frame_color = color::black;
  double coordinate_margin = 0.0;
  double top_margin = 50.0;
  double margin = 50.0;
  double spacing1 = 0.0;
  double inner_frame_width = 0.04;
  double spacing2 = 0.05;
  double outer_frame_width = 0.08;
  Color shading_color = color::gray;
  double shading_line_width = 2.0;
};

class ChessDiagram : public MultiRegion
{
 private:
  ChessDiagramStyle style_;
  cairowindow::Rectangle geometry_;             // The geometry passed to the constructor.
  std::vector<std::shared_ptr<LayerRegion>> regions_;

 public:
  ChessDiagram(cairowindow::Rectangle const& geometry, ChessDiagramStyle style) :
    MultiRegion(style.frame_color, style.shading_line_width), style_(style), geometry_(geometry)
  {
    DoutEntering(dc::notice, "ChessDiagram::ChessDiagram(" << geometry << ", ...)");
  }

  ChessDiagramStyle const& style() const { return style_; }
  cairowindow::Rectangle const& geometry() const { return geometry_; }

 private:
  void draw_regions_on(Layer* layer) override;
};

} // namespace cairowindow::draw

