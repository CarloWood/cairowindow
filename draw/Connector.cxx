#include "sys.h"
#include "Connector.h"

namespace cairowindow::draw {

void Connector::draw_arrow_heads(boost::intrusive_ptr<Layer> const& layer)
{
  if (arrow_head_shape_from_ != cairowindow::Connector::no_arrow)
  {
    layer->draw(arrow_head_from_);
    update_from(arrow_head_from_->arrow_overshoot());
  }
  if (arrow_head_shape_to_ != cairowindow::Connector::no_arrow)
  {
    layer->draw(arrow_head_to_);
    update_to(arrow_head_to_->arrow_overshoot());
  }
}

void Connector::update_from(double overshoot)
{
  double fraction = overshoot / length();
  x1_ += fraction * (x2_ - x1_);
  y1_ += fraction * (y2_ - y1_);
}

void Connector::update_to(double overshoot)
{
  double fraction = overshoot / length();
  x2_ -= fraction * (x2_ - x1_);
  y2_ -= fraction * (y2_ - y1_);
}

} // namespace cairowindow::draw
