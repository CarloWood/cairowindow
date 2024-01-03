#include "sys.h"
#include "Connector.h"

namespace cairowindow::draw {

void Connector::draw_arrow_heads(boost::intrusive_ptr<Layer> const& layer)
{
  if (arrow_head_shape_from_ != cairowindow::Connector::no_arrow)
    layer->draw(arrow_head_from_);
  if (arrow_head_shape_to_ != cairowindow::Connector::no_arrow)
    layer->draw(arrow_head_to_);
}

} // namespace cairowindow::draw
