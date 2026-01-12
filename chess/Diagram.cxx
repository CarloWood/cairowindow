#include "sys.h"
#include "cairowindow/Layer.h"
#include "Diagram.h"
#include "cairowindow/draw/ChessPiece.h"
#include <cctype>

namespace cairowindow::chess {

// Calculate the diagram geometry from the given geometry and style.
//
// The width and height of this geometry is used for the complete widget and therefore
// includes an optional title, a frame around the board and coordinates.
//
// The board exists of a grid of 8x8 squares, which in order to look the same need to
// be at an integer number of pixels distance (on the screen), i.e. the coordinates of
// the returned Rectangle.
//
// All margins specified in the Style are *minimal* margins; they can become
// larger because of rounding off and/or because the geometry passed to this
// function is simply too large in one direction (the board will be centered
// if ChessDiagramStyle::coordinate_margin equals ChessDiagramStyle::margin
// (or not used because no coordinates are drawn) and ChessDiagramStyle::top_margin
// equals ChessDiagramStyle::margin).
//
// However, any excess margin is as much as possible removed by returning
// it as an offset relative to the top-left corner of the specified geometry.
//
// For example:
//
//   Input geometry:
//   +--------------------------------------------------------
//   |
//   |  Output geometry with *extra* offset:
//   |  +---------------------------------
//   |  |                         ^
//   |  |                         | specified top_margin
//   |  |                         |
//   |  |                         v
//   |  |                      +-----------
//   |  |<--specified margin-->| A8 square
//   |  |                      |
//
cairowindow::Geometry Diagram::calculate_geometry(cairowindow::Geometry const& geometry, draw::ChessDiagramStyle const& style)
{
  DoutEntering(dc::notice, "Diagram::calculate_geometry(" << geometry << ", style)");

  bool have_coordinates = style.coordinate_margin() != 0.0;

  double const left_margin = have_coordinates ? style.coordinate_margin() : style.margin();
  double const right_margin = style.margin();
  double const top_margin = style.top_margin();
  double const bottom_margin = have_coordinates ? style.coordinate_margin() : style.margin();
  double const max_horizontal_square_size = (geometry.width() - (left_margin + right_margin)) / 8.0;
  double const max_vertical_square_size = (geometry.height() - (top_margin + bottom_margin)) / 8.0;

  int const square_size_ = std::min(max_vertical_square_size, max_horizontal_square_size);
  int const board_size = 8 * square_size_;

  double output_width = left_margin + board_size + right_margin;
  double output_height = top_margin + board_size + bottom_margin;

  // The center of the input geometry more or less coincides with the center of the output geometry.
  // Calculate the center of the input geometry.
  double geometry_center_x = geometry.offset_x() + 0.5 * geometry.width();
  double geometry_center_y = geometry.offset_y() + 0.5 * geometry.height();
  // Now assume this is the center of the output geometry and find the top-left corner of the chess board.
  // Round this up to the nearest integer.
  int top_left_a8_x = std::ceil(geometry_center_x - 0.5 * output_width + left_margin);
  int top_left_a8_y = std::ceil(geometry_center_y - 0.5 * output_height + top_margin);
  // Store the geometry of the board.
  chess_diagram_.set_board_geometry(top_left_a8_x, top_left_a8_y, board_size);
  // From here, find the top-left corner of the output geometry.
  double top_left_x = top_left_a8_x - left_margin;
  double top_left_y = top_left_a8_y - top_margin;

  Dout(dc::notice, "Returning {" << top_left_x << ", " << top_left_y << ", " << output_width << ", " << output_height << "}");
  return { top_left_x, top_left_y, output_width, output_height };
}

void Diagram::add_to(boost::intrusive_ptr<Layer> const& layer)
{
  // Draw coordinates, frame(s) and shading.
  draw_multi_region_on(layer, &chess_diagram_);

  // Draw the title if any.
  if (title_)
  {
    title_->move_to(chess_diagram_.geometry().offset_x() + 0.5 * chess_diagram_.geometry().width(),
        chess_diagram_.geometry().offset_y() + 0.5 * (chess_diagram_.style().top_margin() -
          chess_diagram_.frame_thickness()) - title_->style().offset());
    draw_layer_region_on(layer, title_);
  }

  // Register this diagram and its geometry with the associated Window so that we can find which printable is under the mouse if needed.
  layer->window()->add_printable(this);
}

void Diagram::add_piece(boost::intrusive_ptr<Layer> const& layer,
    draw::ChessPieceStyle const& style,
    chess::Piece const& piece)
{
  int col = piece.col();
  int row = piece.row();

  // Calculate the geometry of the square that this piece is placed on.
  int const square_size = chess_diagram_.board_size() / 8;
  double const top_left_x = chess_diagram_.top_left_x() + col * square_size;
  double const top_left_y = chess_diagram_.top_left_y() + (7 - row) * square_size;

  piece.draw_object_ = std::make_shared<draw::ChessPiece>(
      top_left_x, top_left_y, square_size,
      piece.color(), piece.piece(), style);
  draw_layer_region_on(layer, piece.draw_object_);
}

void Diagram::place_piece(boost::intrusive_ptr<Layer> const& layer,
    EColor color, EPiece piece, int col, int row, draw::ChessPieceStyle const& style)
{
  pieces_.emplace_back(color, piece, col, row);
  add_piece(layer, style, pieces_.back());
}

void Diagram::clear()
{
  pieces_.clear();
}

bool Diagram::load_FEN(boost::intrusive_ptr<Layer> const& layer, std::string const& FEN, draw::ChessPieceStyle const& style)
{
  DoutEntering(dc::notice, "Diagram::load_FEN(\"" << FEN << "\")");

  // Clear the position.
  clear();
  std::string::const_iterator iter = FEN.begin();
  char c;
  // Eat preceding spaces.
  while (*iter == ' ')
    ++iter;
  if (iter == FEN.end())
    return false;
  // Field 1: Piece placement.
  EColor color;
  int col = 0, row = 7;
  while ((c = *iter++) != ' ')
  {
    color = white;
    switch (c)
    {
      case '/':
        if (col != 8 || row <= 0)
	  return false;
        col = 0;
        --row;
        break;
      case '1':
      case '2':
      case '3':
      case '4':
      case '5':
      case '6':
      case '7':
      case '8':
        col += (int)(c - '0');
	break;
      case 'p':
      case 'r':
      case 'n':
      case 'b':
      case 'q':
      case 'k':
        color = black;
        /* FALL-THROUGH */
      case 'P':
      case 'R':
      case 'N':
      case 'B':
      case 'Q':
      case 'K':
      {
        c = std::toupper(c);
	EPiece type;
	if (c == 'P')
	  type = pawn;
	else if (c == 'R')
	  type = rook;
	else if (c == 'N')
	  type = knight;
	else if (c == 'B')
	  type = bishop;
	else if (c == 'Q')
	  type = queen;
	else // if (c == 'K')
	  type = king;
	place_piece(layer, color, type, col, row, style);
        ++col;
        break;
      }
      default:
        return false;
    }
    if (iter == FEN.end())
      return false;
  }
  if (col != 8)
    return false;
  if (iter == FEN.end())
    return false;
  // Eat extra spaces.
  while (*iter == ' ')
    if (++iter == FEN.end())
      return false;
  // Field 2: Active color.
  c = *iter++;
  if (c != 'w' && c != 'b')
    return false;
  to_move_ = (c == 'w') ? white : black;
  if (iter == FEN.end() || *iter != ' ')
    return false;
  // Eat space and possibly extra spaces.
  while (*iter == ' ')
    if (++iter == FEN.end())
      return false;
#if 0
  // Field 3: Castling availability.
  uint8_t white_castle_flags = white_rook_queen_side_moved | white_king_moved | white_rook_king_side_moved;
  uint8_t black_castle_flags = black_rook_queen_side_moved | black_king_moved | black_rook_king_side_moved;
#endif
  while ((c = *iter++) != ' ')
  {
    if (c == '-')
      break;
    switch (c)
    {
      case 'K':
//        white_castle_flags &= ~(white_king_moved | white_rook_king_side_moved);
	break;
      case 'Q':
//        white_castle_flags &= ~(white_king_moved | white_rook_queen_side_moved);
	break;
      case 'k':
//        black_castle_flags &= ~(black_king_moved | black_rook_king_side_moved);
	break;
      case 'q':
//        black_castle_flags &= ~(black_king_moved | black_rook_queen_side_moved);
	break;
      default:
        return false;
    }
    if (iter == FEN.end())
      return false;
  }
//  castle_flags = white_castle_flags | black_castle_flags;
  if (iter == FEN.end())
    return false;
  // Eat extra spaces.
  while (*iter == ' ')
    if (++iter == FEN.end())
      return false;
  // Field 4: En passant target square in algebraic notation.
  if ((c = *iter++) != '-')
  {
    if (c < 'a' || c > 'h' || iter == FEN.end())
      return false;
    int col = c - 'a';
    c = *iter++;
    if (c < '1' || c > '8' || iter == FEN.end())
      return false;
    int row = c - '1';
#if 0
    Code other_pawn;
    if (M_to_move == white)
    {
      other_pawn = white_pawn;
      if (row != 5 || M_pieces[Index(col, row - 1)] != black_pawn || M_pieces[Index(col, row)] != nothing)
        return false;
    }
    else
    {
      other_pawn = black_pawn;
      if (row != 2 || M_pieces[Index(col, row + 1)] != white_pawn || M_pieces[Index(col, row)] != nothing)
        return false;
    }
    set_en_passant(Index(col, row));
#endif
  }
  if (iter == FEN.end() || *iter != ' ')
    return false;
  // Eat space and possibly extra spaces.
  while (*iter == ' ')
    if (++iter == FEN.end())
      return false;
  // Field 5: Halfmove clock: This is the number of halfmoves since the last pawn advance or capture.
//  M_half_move_clock = 0;
  while ((c = *iter++) != ' ')
  {
    if (!std::isdigit(c))
      return false;
//    M_half_move_clock = 10 * M_half_move_clock + c - '0';
    if (iter == FEN.end())
      return false;
  }
  // Eat extra spaces.
  while (*iter == ' ')
    if (++iter == FEN.end())
      return false;
  // Field 6: Fullmove number: The number of the full move.
  int M_full_move_number = 0;
  while (iter != FEN.end() && *iter != ' ')
  {
    if (!std::isdigit(*iter))
      return false;
    M_full_move_number = 10 * M_full_move_number + *iter++ - '0';
  }
  if (M_full_move_number == 0)
    return false;
  // Success.
  return true;
}

} // namespace cairowindow::chess
