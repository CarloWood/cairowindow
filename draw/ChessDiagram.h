#pragma once

#include "MultiRegion.h"
#include "Line.h"

namespace cairowindow::draw {

#define cairowindow_ChessDiagram_FOREACH_MEMBER(X, ...) \
  X(Color, frame_color, Color{}, __VA_ARGS__) \
  X(double, coordinate_margin, -1.0, __VA_ARGS__) \
  X(double, top_margin, -1.0, __VA_ARGS__) \
  X(double, margin, -1.0, __VA_ARGS__) \
  X(double, spacing1, -1.0, __VA_ARGS__) \
  X(double, inner_frame_width, -1.0, __VA_ARGS__) \
  X(double, spacing2, -1.0, __VA_ARGS__) \
  X(double, outer_frame_width, -1.0, __VA_ARGS__) \
  X(Color, shading_color, Color{}, __VA_ARGS__) \
  X(double, shading_line_width, -1.0, __VA_ARGS__)

#define cairowindow_ChessDiagram_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_ChessDiagram_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for ChessDiagramStyle.
struct ChessDiagramStyleParamsDefault
{
  static constexpr Color frame_color = color::black;
  static constexpr double coordinate_margin = 0.0;
  static constexpr double top_margin = 50.0;
  static constexpr double margin = 50.0;
  static constexpr double spacing1 = 0.0;
  static constexpr double inner_frame_width = 0.04;
  static constexpr double spacing2 = 0.05;
  static constexpr double outer_frame_width = 0.08;
  static constexpr Color shading_color = color::gray;
  static constexpr double shading_line_width = 2.0;
};

// Declare ChessDiagramStyle.
DECLARE_STYLE(ChessDiagram, ChessDiagramStyleParamsDefault);

class ChessDiagram : public MultiRegion
{
 private:
  ChessDiagramStyle style_;
  cairowindow::Geometry geometry_;             // The geometry passed to the constructor.

  // These are calculated by the call to Diagram::add_to.
  int top_left_x_;                              // x coordinate of top-left corner of the a8 square.
  int top_left_y_;                              // y coordinate of top-left corner of the a8 square.
  int board_size_;                              // The width of the chess board (this is divisible by 8).
  double frame_thickness_{};
  int square_size_{};
 std::vector<std::shared_ptr<LayerRegion>> regions_;   // The LayerRegion's used, filled by draw_regions_on.

 public:
  ChessDiagram(cairowindow::Geometry const& geometry, ChessDiagramStyle const& style) :
    MultiRegion(style.frame_color(), style.shading_line_width()), style_(style), geometry_(geometry)
  {
    DoutEntering(dc::notice, "ChessDiagram::ChessDiagram(" << geometry << ", ...)");
  }

  ChessDiagramStyle const& style() const { return style_; }
  cairowindow::Geometry const& geometry() const { return geometry_; }

  double frame_thickness() const { ASSERT(frame_thickness_ > 0.0); return frame_thickness_; }
  int square_size() const { ASSERT(square_size_ > 0); return square_size_; }

  void set_board_geometry(int top_left_x, int top_left_y, int board_size)
  {
    ASSERT(board_size % 8 == 0);
    top_left_x_ = top_left_x;
    top_left_y_ = top_left_y;
    board_size_ = board_size;
  }

  int top_left_x() const { return top_left_x_; }
  int top_left_y() const { return top_left_y_; }
  int board_size() const { return board_size_; }

 private:
  void draw_regions_on(Layer* layer) override;
};

} // namespace cairowindow::draw
