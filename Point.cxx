#include "sys.h"
#include "Plot.h"
#include "Layer.h"

namespace cairowindow::plot::cs {

// Definition of the specialization for cs = CS::plot.

template<>
void Point<CS::plot>::moved(Plot* plot, cairowindow::cs::Point<CS::plot> const& new_position)
{
  *this = new_position;
  plot->add_point(draw_object_->layer(), draw_object_->point_style(), *this);
}

template<>
void Point<CS::plot>::move(Plot& UNUSED_ARG(plot), cairowindow::cs::Point<CS::plot> const& new_position)
{
  Layer* layer = draw_object_->layer();
  Window* window = layer->window();
  window->move_draggable(this, index_, new_position);
}

} // namespace cairowindow::plot::cs
