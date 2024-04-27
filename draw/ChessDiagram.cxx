#include "sys.h"
#include "ChessDiagram.h"
#include "Rectangle.h"
#include "Text.h"
#include "cairowindow/Layer.h"
#ifdef CWDEBUG
#include "cairowindow/debugcairo.h"
#endif

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
// |       |<-|-frame_thickness---->|
// |<---------|-left_margin|------> |                                   ChessDiagramStyle::coordinate_margin (if coordinates must be drawn).
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
  int const square_size = board_size_ >> 3;
  int const bottom_right_x = top_left_x_ + board_size_;
  int const bottom_right_y = top_left_y_ + board_size_;

  for (int frame = 0; frame <= (style_.outer_frame_width() == 0.0 ? 0 : 1); ++frame)
  {
    double offset = (style_.spacing1() +
      ((frame == 0) ? 0.5 * style_.inner_frame_width() : style_.inner_frame_width() + style_.spacing2() + 0.5 * style_.outer_frame_width())) *
      square_size;
    double x1 = top_left_x_ - 0.5 - offset;
    double y1 = top_left_y_ - 0.5 - offset;
    double x2 = bottom_right_x + 0.5 + offset;
    double y2 = bottom_right_y + 0.5 + offset;
    double frame_width = ((frame == 0) ? style_.inner_frame_width() : style_.outer_frame_width()) * square_size;
    LineCap line_cap = (frame == 0) ? LineCap::square : LineCap::round;
    regions_.emplace_back(std::make_shared<Line>(x1, y1 + 0.2, x1, y2 - 0.2,
          LineStyle{{.line_color = color_, .line_width = frame_width, .line_cap = line_cap}}));
    regions_.emplace_back(std::make_shared<Line>(x2, y1 + 0.2, x2, y2 - 0.2,
          LineStyle{{.line_color = color_, .line_width = frame_width, .line_cap = line_cap}}));
    regions_.emplace_back(std::make_shared<Line>(x1, y1, x2, y1,
          LineStyle{{.line_color = color_, .line_width = frame_width}}));
    regions_.emplace_back(std::make_shared<Line>(x1, y2, x2, y2,
          LineStyle{{.line_color = color_, .line_width = frame_width}}));
  }

  // Shading.
  double y1 = bottom_right_y;
  for (int yi = 0; yi < 8; ++yi)
  {
    double x1 = top_left_x_ + ((yi & 1) ? square_size : 0.0);
    for (int xi = yi & 1; xi < 8; xi += 2)
    {
      double x2 = x1 + square_size;
      double y2 = y1 - square_size;
      LineStyle shading_line_style({.line_color = style_.shading_color(), .line_width = 0.01 * square_size * line_width_});
      for (double sh = 0.0; sh < square_size; sh += square_size / 7.5)
      {
        regions_.emplace_back(std::make_shared<Line>(x1 + sh, y1, x2, y2 + sh, shading_line_style));
        if (sh > 0.0)
          regions_.emplace_back(std::make_shared<Line>(x1, y1 - sh, x2 - sh, y2, shading_line_style));
      }
      x1 += 2 * square_size;
    }
    y1 -= square_size;
  }

  frame_thickness_ = (style_.spacing1() + style_.inner_frame_width() +
    ((style_.outer_frame_width() == 0.0) ? 0.0 : style_.spacing2() + style_.outer_frame_width())) * square_size;

  if (style_.coordinate_margin() != 0.0)
  {
#ifdef CWDEBUG
    using namespace debugcairo;
#endif
    double x1 = top_left_x_;
    double y1 = top_left_y_ + 0.5 * square_size;
    double offset = frame_thickness_ + 0.15 * square_size;
    TextStyle coordinates_style({.position = centered_left_of, .font_size = 0.5 * square_size, .offset = offset});
    for (int yi = 0; yi < 8; ++yi)
    {
      regions_.emplace_back(std::make_shared<Text>(std::to_string(8 - yi), x1, y1, coordinates_style));
      y1 += square_size;
    }
    x1 += 0.5 * square_size;
    y1 = top_left_y_ + board_size_;
    coordinates_style.setup(layer->cr());
    cairo_text_extents_t extents;
    cairo_text_extents(layer->cr(), "abcdefgh", &extents);
    double y_bearing = extents.y_bearing;
    TextStyle a_h_style = coordinates_style({.position = centered_below_no_bearing, .offset = offset - y_bearing});
    for (char xi = 'a'; xi <= 'h'; ++xi)
    {
      std::string coordinate(1, xi);
      regions_.emplace_back(std::make_shared<Text>(coordinate, x1, y1, a_h_style));
      x1 += square_size;
    }
  }

  for (auto const& region : regions_)
    layer->draw(region);
}

} // namespace cairowindow::draw
