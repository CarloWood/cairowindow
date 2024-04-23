#include "sys.h"
#include "cairowindow/Layer.h"
#include "Diagram.h"

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
cairowindow::Rectangle Diagram::calculate_geometry(cairowindow::Rectangle const& geometry, draw::ChessDiagramStyle const& style)
{
  DoutEntering(dc::notice, "Diagram::calculate_geometry(" << geometry << ", style)");

  bool have_coordinates = style.coordinate_margin != 0.0;

  double const left_margin = have_coordinates ? style.coordinate_margin : style.margin;
  double const right_margin = style.margin;
  double const top_margin = style.top_margin;
  double const bottom_margin = style.margin;
  double const max_horizontal_square_size = (geometry.width() - (left_margin + right_margin)) / 8.0;
  double const max_vertical_square_size = (geometry.height() - (top_margin + bottom_margin)) / 8.0;

  int const square_size = std::min(max_vertical_square_size, max_horizontal_square_size);
  int const board_size = 8 * square_size;

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
  // From here, find the top-left corner of the output geometry.
  double top_left_x = top_left_a8_x - left_margin;
  double top_left_y = top_left_a8_y - top_margin;

  Dout(dc::notice, "Returning {" << top_left_x << ", " << top_left_y << ", " << output_width << ", " << output_height << "}");
  return { top_left_x, top_left_y, output_width, output_height };
}

void Diagram::add_to(boost::intrusive_ptr<Layer> const& layer)
{
  // Set a title.
  if (title_)
  {
    title_->move_to(chess_diagram_.geometry().offset_x() + 0.5 * chess_diagram_.geometry().width(),
        chess_diagram_.geometry().offset_y() + 0.5 * chess_diagram_.style().top_margin - title_->style().offset());
    draw_layer_region_on(layer, title_);
  }

  // Set ranges on the plot area and draw it.
  draw_multi_region_on(layer, &chess_diagram_);

  // Draw coordinates.
  //...

  // Register this diagram and its geometry with the associated Window so that we can find which printable is under the mouse if needed.
  layer->window()->add_printable(this);
}

} // namespace cairowindow::chess
