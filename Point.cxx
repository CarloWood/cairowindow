#include "sys.h"
#include "Point.h"
#include "Direction.h"
#include "Vector.h"

namespace cairowindow {

Point Point::operator+(Direction const& direction)
{
  return {x_ + direction.x(), y_ + direction.y()};
}

Point Point::operator+(Vector const& v)
{
  return {x_ + v.x(), y_ + v.y()};
}

Point Point::operator-(Vector const& v)
{
  return {x_ - v.x(), y_ - v.y()};
}

Point& Point::operator+=(Direction const& direction)
{
  x_ += direction.x();
  y_ += direction.y();
  return *this;
}

Point& Point::operator+=(Vector const& v)
{
  x_ += v.x();
  y_ += v.y();
  return *this;
}

Point& Point::operator-=(Vector const& v)
{
  x_ -= v.x();
  y_ -= v.y();
  return *this;
}

Vector operator-(Point const& from, Point const& to)
{
  return {from, to};
}

#ifdef CWDEBUG
void Point::print_on(std::ostream& os) const
{
  os << "{x_:" << x_ << ", y_:" << y_ << '}';
}
#endif

} // namespace cairowindow
