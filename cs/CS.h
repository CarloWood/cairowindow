#pragma once

namespace cairowindow {

enum class CS
{
  painter,              // The space defined by the painter’s Current Transformation Matrix (CTM) at the instant you issue a painter->draw*() call.
  centered,             // Coordinate system with origin in the middle of the window, where -1 corresponds with the bottom of the window and 1 with the top.
  pixels,               // The final coordinate system of the window in pixels.
  plot
};

} // namespace cairowindow
