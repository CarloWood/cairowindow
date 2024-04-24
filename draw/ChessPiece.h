#pragma once

#include "cairowindow/LayerRegion.h"
#include "cairowindow/chess/Diagram.h"  // For EPiece and EColor.

namespace cairowindow::draw {

#define cairowindow_ChessPiece_FOREACH_MEMBER(X, ...) \
  X(Color, line_color, Color{}, __VA_ARGS__)

#define cairowindow_ChessPiece_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_ChessPiece_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for ChessPieceStyle.
struct ChessPieceStyleParamsDefault
{
  static constexpr Color line_color = color::black;
};

// Declare ChessPieceStyle.
DECLARE_STYLE(ChessPiece, ChessPieceStyleParamsDefault);

class ChessPiece : public LayerRegion
{
 private:
  double x1_;
  double y1_;
  int square_size_;             // The size of the piece.
  chess::EColor color_;
  chess::EPiece piece_;
  ChessPieceStyle style_;

 public:
  ChessPiece(double x1, double y1, int square_size, chess::EColor color, chess::EPiece piece, ChessPieceStyle const& style) :
    x1_(x1), y1_(y1), square_size_(square_size), color_(color), piece_(piece), style_(style) { }

 private:
  StrokeExtents do_draw(cairo_t* cr) override;
};

} // namespace cairowindow::draw
