#include "sys.h"
#include "SampleNode.h"

namespace gradient_descent {

#ifdef CWDEBUG
void LocalExtreme::debug_print_label(char const* left_or_right, const_iterator neighbor) const
{
  Dout(dc::notice, "Setting " << left_or_right << " neighbor of \"" << label() << "] to [" << neighbor->local_extreme().label() << "\".");
}
#endif

} // namespace gradient_descent
