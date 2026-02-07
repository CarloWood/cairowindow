#pragma once

#include "Rectangle.h"
#include "math/cs/TransformOperators.h"

namespace cairowindow::cs {

template<CS from_cs, CS to_cs, bool inverted, math::AffineTransformConcept AffineTransformBackend>
Rectangle<to_cs> operator*(Rectangle<from_cs> const& rectangle, math::Transform<from_cs, to_cs, inverted, AffineTransformBackend> const& transform)
{
  math::cs::Point<from_cs> p_from_cs{rectangle.offset_x(), rectangle.offset_y()};
  math::cs::Size<from_cs> s_from_cs{rectangle.width(), rectangle.height()};

  math::cs::Point<to_cs> p_to_cs = p_from_cs * transform;
  math::cs::Size<to_cs> s_to_cs = s_from_cs * transform;

  return {p_to_cs, s_to_cs};
}

} // namespace cairowindow::cs
