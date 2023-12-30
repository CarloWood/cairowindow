#pragma once

#include "cairowindow/MultiRegion.h"
#include "Line.h"
#include <memory>
#include <vector>

namespace cairowindow::draw {

class Curve : public MultiRegion
{
 private:
  std::vector<std::shared_ptr<Line>> lines_;

  // Implementation of MultiRegion.
  void draw_regions_on(Layer* layer) override;

 public:
  Curve(LineStyle line_style) :
    MultiRegion(line_style.line_color, line_style.line_width) { }

  Curve(std::vector<std::shared_ptr<Line>> const& lines, LineStyle line_style) :
    MultiRegion(line_style.line_color, line_style.line_width), lines_(lines) { }

  std::vector<std::shared_ptr<Line>> const& lines() const { return lines_; }
  std::vector<std::shared_ptr<Line>>& lines() { return lines_; }
};

} // namespace cairowindow::draw
