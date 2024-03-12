#pragma once

#include "Draggable.h"
#include <memory>
#ifdef CWDEBUG
#include "utils/has_print_on.h"
#endif

namespace cairowindow {
#ifdef CWDEBUG
using utils::has_print_on::operator<<;
#endif

class Direction;
class Vector;

class Point
{
 private:
  double x_;
  double y_;

 public:
  Point() = default;
  Point(double x, double y) : x_(x), y_(y) { }

  double x() const { return x_; }
  double y() const { return y_; }

  Point operator+(Direction const& direction);
  Point operator+(Vector const& v);
  Point operator-(Vector const& v);
  Point& operator+=(Direction const& direction);
  Point& operator+=(Vector const& v);
  Point& operator-=(Vector const& v);
  friend Vector operator-(Point const& to, Point const& from);
  friend bool operator!=(Point const& p1, Point const& p2);

#ifdef CWDEBUG
  void print_on(std::ostream& os) const;
#endif
};

namespace draw {
class Point;
} // namespace draw

namespace plot {
class Plot;

class Point : public cairowindow::Point, public Draggable
{
 public:
  Point() = default;
  Point(cairowindow::Point const& point) : cairowindow::Point(point) { }

 private:
  friend class Plot;
  mutable std::shared_ptr<draw::Point> draw_object_;

  // Implementation of Draggable.
  cairowindow::Rectangle const& geometry() const override;
  void moved(Plot* plot, cairowindow::Point const& new_position) override;

#ifdef CWDEBUG
 public:
  void print_on(std::ostream& os) const override { cairowindow::Point::print_on(os); }
#endif
};

} // namespace plot
} // namespace cairowindow
