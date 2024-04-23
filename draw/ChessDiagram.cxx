#include "sys.h"
#include "ChessDiagram.h"
#include "Rectangle.h"
#include "cairowindow/Layer.h"

namespace cairowindow::draw {

// +-----------------------------------------------------------------------------
// |                                                                           ^
// |                 Title here, for which we need a vertical margin           |
// |       +-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,____ outer_frame_width   | ChessDiagramStyle::top_margin
// |       |X+'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-                         |
// |       |X| +-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,-,____ inner_frame_width   |
// |       |X| |X+'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-'-                         v
// |       |X| |X| +------+------+------+------+------
// |       |X| |X| |      |      |      |      |
// |    8  |X| |X| |      |      |      |      |
// |       |X| |X| |      |      |      |      |
// |       |X| |X| +------+------+------+------+------
// |    7  |X| |X| |      |      |      |      |
// |       |X| |X| |      |      |      |      |
// |    .  |X| |X| |      |      |      |      |
// |    .  |X| |X| +------+------+------+------+------
// |<----->|X| |X| |      |      |      |      |
// | |              <----->
// |  \_ left_margin  \_ square_size (must be an integer)
//
// Zooming in on the pixels one side of the board:
//
//              __ outer_frame_width (3 pixels in this case).                   ChessDiagramStyle::outer_frame_width (or 0 when no frame).
//             /            __ inner_frame_width (2 pixels in this case).       ChessDiagramStyle::inner_frame_width.
//            /            /
//         <----->       <--->     _spacing1 (3 pixels here)                    ChessDiagramStyle::spacing1
//                    _spacing2   /                                             ChessDiagramStyle::spacing2
//                   / (4 pixels)/                                                      ^
//         ├─┼─┼─┤  /    ├─┼─┤  /  ├─┼─                                                 |
//         ├─┼─┼─┤<----->├─┼─┤<--->├─┼─                                             These style values are
//         ├─┼─┼─┤       ├─┼─┤     ├─┼─ pixels used                                 specified as fraction of
//         ├─┼─┼─┤       ├─┼─┤     ├─┼─ by square in                                square_size.
//         ├─┼─┼─┤       ├─┼─┤     ├─┼─ file a.
//         ├─┼─┼─┤       ├─┼─┤     ├─┼─
//         ├─┼─┼─┤       ├─┼─┤     ├─┼─
//         ├─┼─┼─┤       ├─┼─┤     ├─┼─
//            ^            ^        ^
// |<---------|-left-margin|------> |                                   ChessDiagramStyle::coordinate_margin (if coordinates must be drawn).
//            |            |        |                                   ChessDiagramStyle::margin (when there are no coordinates).
//            |            |<-----> |   - inner_frame_distance
//            |<-----------|------> |   - outer_frame_distance
//            |            |        |
//        x-coordinate x-coordinate x-coordinate of first pixel
//        outer frame  inner frame  of square (51, so that the
//        line         line         left-margin here is 50.5).
//        (e.g 40)     (46.5)       This must be an integer.
//
void ChessDiagram::draw_regions_on(Layer* layer)
{
  constexpr int horizontal = 0;
  constexpr int vertical = 1;

  bool have_coordinates = style_.coordinate_margin != 0.0;
  double const left_margin = have_coordinates ? style_.coordinate_margin : style_.margin;
  double const right_margin = style_.margin;
  double const top_margin = style_.top_margin;
  double const bottom_margin = style_.margin;

  int const bottom_left_x = std::round(geometry_.offset_x() + right_margin);
  int const bottom_left_y = std::round(geometry_.offset_y() + geometry_.height() - bottom_margin);
  int const board_size = std::round((geometry_.width() - left_margin - right_margin));
  ASSERT(board_size % 8 == 0);
  int const square_size = board_size >> 3;

  int const top_right_x = bottom_left_x + board_size;
  int const top_right_y = bottom_left_y - board_size;

  for (int frame = 0; frame <= (style_.outer_frame_width == 0.0 ? 0 : 1); ++frame)
  {
    double offset = (style_.spacing1 +
      ((frame == 0) ? 0.5 * style_.inner_frame_width : style_.inner_frame_width + style_.spacing2 + 0.5 * style_.outer_frame_width)) *
      square_size;
    double x1 = bottom_left_x - 0.5 - offset;
    double y1 = bottom_left_y + 0.5 + offset;
    double x2 = top_right_x + 0.5 + offset;
    double y2 = top_right_y - 0.5 - offset;
    double frame_width = ((frame == 0) ? style_.inner_frame_width : style_.outer_frame_width) * square_size;
    LineCap line_cap = (frame == 0) ? LineCap::square : LineCap::round;
    regions_.emplace_back(std::make_shared<Line>(x1, y2 + 0.2, x1, y1 - 0.2,
          LineStyle{{.line_color = color_, .line_width = frame_width, .line_cap = line_cap}}));
    regions_.emplace_back(std::make_shared<Line>(x2, y2 + 0.2, x2, y1 - 0.2,
          LineStyle{{.line_color = color_, .line_width = frame_width, .line_cap = line_cap}}));
    regions_.emplace_back(std::make_shared<Line>(x1, y2, x2, y2,
          LineStyle{{.line_color = color_, .line_width = frame_width}}));
    regions_.emplace_back(std::make_shared<Line>(x1, y1, x2, y1,
          LineStyle{{.line_color = color_, .line_width = frame_width}}));
  }

  double y1 = bottom_left_y;
  for (int yi = 0; yi < 8; ++yi)
  {
    double x1 = bottom_left_x + ((yi & 1) ? square_size : 0.0);
    for (int xi = yi & 1; xi < 8; xi += 2)
    {
      double x2 = x1 + square_size;
      double y2 = y1 - square_size;
      for (double sh = 0.0; sh < square_size; sh += square_size / 7.5)
      {
        regions_.emplace_back(std::make_shared<Line>(x1 + sh, y1, x2, y2 + sh,
              LineStyle{{.line_color = style_.shading_color, .line_width = line_width_}}));
        if (sh > 0.0)
          regions_.emplace_back(std::make_shared<Line>(x1, y1 - sh, x2 - sh, y2,
                LineStyle{{.line_color = style_.shading_color, .line_width = line_width_}}));
      }
      x1 += 2 * square_size;
    }
    y1 -= square_size;
  }

  for (auto const& region : regions_)
    layer->draw(region);
}

} // namespace cairowindow::draw
