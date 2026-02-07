#pragma once

#include "Point.h"
#include "Size.h"
#include "math/Transform.h"

namespace cairowindow::cs {

template<CS from_cs, CS to_cs, bool inverted, math::AffineTransformConcept AffineTransformBackend>
Point<to_cs> operator*(Point<from_cs> const& point, math::Transform<from_cs, to_cs, inverted, AffineTransformBackend> const& transform)
{
  if constexpr (!inverted)
  {
    auto const [x, y] = transform.map_point(point.x(), point.y());
    return {x, y};
  }
  else
  {
    auto const inv = transform.inverted();
    auto const [x, y] = inv.map_point(point.x(), point.y());
    return {x, y};
  }
}

template<CS from_cs, CS to_cs, bool inverted, math::AffineTransformConcept AffineTransformBackend>
Size<to_cs> operator*(Size<from_cs> const& size, math::Transform<from_cs, to_cs, inverted, AffineTransformBackend> const& transform)
{
  // Just scale; scale the X and Y axis vectors by the full linear part.
  if constexpr (!inverted)
  {
    auto const [sx, sy] = transform.scale_factors();
    return {size.width() * sx, size.height() * sy};
  }
  else
  {
    auto const inv = transform.inverted();
    auto const [sx, sy] = inv.scale_factors();
    return {size.width() * sx, size.height() * sy};
  }
}

#if 0 //FIXME: move this to cairowindow
template<CS from_cs, CS to_cs, bool inverted, math::AffineTransformConcept AffineTransformBackend>
Rectangle<to_cs> operator*(Rectangle<from_cs> const& rectangle, math::Transform<from_cs, to_cs, inverted, AffineTransformBackend> const& transform)
{
  Point<from_cs> p_from_cs{rectangle.offset_x(), rectangle.offset_y()};
  Size<from_cs> s_from_cs{rectangle.width(), rectangle.height()};

  Point<to_cs> p_to_cs = p_from_cs * transform;
  Size<to_cs> s_to_cs = s_from_cs * transform;

  return {p_to_cs, s_to_cs};
}
#endif

} // namespace cairowindow::cs
