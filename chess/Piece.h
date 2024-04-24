#pragma once

#include "enums.h"
#include <memory>

namespace cairowindow {

// A chess piece, including the chess coordinates where it stands on a board.
// The board itself is not stored, nor does this class contain any information
// on how to draw the piece.
class ChessPiece
{
 private:
  chess::EColor color_;
  chess::EPiece piece_{};
  int col_;                     // 0..7 corresponding with the file A through H.
  int row_;                     // 0..7 corresponding with the row 1 through 8.

 public:
  ChessPiece() = default;
  ChessPiece(chess::EColor color, chess::EPiece piece, int col, int row) :
    color_(color), piece_(piece), col_(col), row_(row) { }

  chess::EPiece piece() const { return piece_; }
  chess::EColor color() const { return color_; }
  int col() const { return col_; }
  int row() const { return row_; }
};

namespace draw {
class ChessPiece;
} // namespace draw

namespace chess {

// Like the above, but including a draw_object_.
// Only a chess::Diagram can create this object; call Diagram::place_piece.
class Piece : public cairowindow::ChessPiece
{
 public:
  using cairowindow::ChessPiece::ChessPiece;
  Piece(cairowindow::ChessPiece const& chess_piece) : cairowindow::ChessPiece(chess_piece) { }

 private:
  friend class Diagram;
  mutable std::shared_ptr<draw::ChessPiece> draw_object_;
};

} // namespace chess
} // namespace cairowindow
