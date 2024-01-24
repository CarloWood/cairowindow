#include "sys.h"
#include "Pixel.h"
#ifdef CWDEBUG
#include <iostream>
#endif

namespace cairowindow {

#ifdef CWDEBUG
void Pixel::print_on(std::ostream& os) const
{
  os << "{x_:" << x_ << ", y_:" << y_ << '}';
}
#endif

} // namespace cairowindow
