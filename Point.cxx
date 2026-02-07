#include "sys.h"
#include "Layer.h"
#include "plot/Point.h"

namespace cairowindow::plot::cs {

// Definition of the specialization for cs = csid::plot.

template<>
void Point<csid::plot>::moved(cairowindow::Point const& new_position)
{
  *this = new_position;
}

template<>
void Point<csid::plot>::move(Plot& UNUSED_ARG(plot), math::cs::Point<csid::plot> const& new_position)
{
  Layer* layer = draw_object_->layer();
  Window* window = layer->window();
  window->move_draggable(this, index_, new_position);
}

} // namespace cairowindow::plot::cs
