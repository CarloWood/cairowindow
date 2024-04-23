#pragma once

#include "cairowindow/Printable.h"
#include "cairowindow/draw/ChessDiagram.h"
#include "cairowindow/Rectangle.h"
#include "cairowindow/draw/Text.h"
#include <boost/intrusive_ptr.hpp>
#include "debug.h"

namespace cairowindow::chess {

struct TitleStyleDefaults : draw::TextStyleParamsDefault
{
  static constexpr draw::TextPosition position = draw::centered;
  static constexpr double font_size = 24.0;
};

} // namespace cairowindow::chess
namespace cairowindow::draw {

#define cairowindow_ChessTitle_FOREACH_MEMBER cairowindow_TextBase_FOREACH_MEMBER
#define cairowindow_ChessTitle_FOREACH_STYLE_MEMBER cairowindow_TextBase_FOREACH_MEMBER

DECLARE_STYLE_WITH_BASE(ChessTitle, TextBase, chess::TitleStyleDefaults)

#undef cairowindow_ChessTitle_FOREACH_MEMBER
#undef cairowindow_ChessTitle_FOREACH_STYLE_MEMBER

} // namespace cairowindow::draw
namespace cairowindow::chess {

class Diagram : public Printable
{
 private:
  draw::ChessDiagram chess_diagram_;
  std::shared_ptr<draw::Text> title_;

 public:
  Diagram(Rectangle const& geometry, draw::ChessDiagramStyle chess_diagram_style, std::string title, draw::ChessTitleStyle title_style) :
    chess_diagram_(calculate_geometry(geometry, chess_diagram_style), chess_diagram_style),
    title_(std::make_shared<draw::Text>(title, chess_diagram_.geometry().offset_x() + 0.5 * chess_diagram_.geometry().width(),
        chess_diagram_.geometry().offset_y() + 0.5 * chess_diagram_style.top_margin() - title_style.offset(), title_style)) { }

  Diagram(Rectangle const& geometry, draw::ChessDiagramStyle chess_diagram_style) :
    chess_diagram_(calculate_geometry(geometry, chess_diagram_style), chess_diagram_style) { }

  Rectangle const& geometry() const override { return chess_diagram_.geometry(); }
  void add_to(boost::intrusive_ptr<Layer> const& layer);

 private:
  Rectangle calculate_geometry(Rectangle const& geometry, draw::ChessDiagramStyle const& style);
};

} // namespace cairowindow::chess
