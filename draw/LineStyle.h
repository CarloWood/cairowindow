#pragma once

#include "cairowindow/Color.h"
#include "cairowindow/Style.h"

namespace cairowindow::draw {

enum class LineCap
{
  undefined, butt, round, square
};

std::ostream& operator<<(std::ostream& os, LineCap line_cap);

#define cairowindow_Line_FOREACH_MEMBER(X, ...) \
  X(Color, line_color, Color{}, __VA_ARGS__) \
  X(double, line_width, -1.0, __VA_ARGS__) \
  X(std::vector<double>, dashes, std::vector<double>{-1.0}, __VA_ARGS__) \
  X(double, dashes_offset, 12345678.9, __VA_ARGS__) \
  X(LineCap, line_cap, LineCap::undefined, __VA_ARGS__)

#define cairowindow_Line_FOREACH_STYLE_MEMBER(X, ...) \
  cairowindow_Line_FOREACH_MEMBER(X, __VA_ARGS__)

// Define default values for LineStyle.
struct LineStyleParamsDefault
{
  static constexpr Color line_color = color::indigo;
  static constexpr double line_width = 2.0;
  static constexpr std::vector<double> dashes = {};
  static constexpr double dashes_offset = 0.0;
  static constexpr LineCap line_cap = LineCap::butt;
};

// Declare LineStyle.
DECLARE_STYLE(Line, LineStyleParamsDefault);

} // namespace cairowindow::draw
