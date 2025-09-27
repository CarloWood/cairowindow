#include "sys.h"
#include "SampleNode.h"

namespace gradient_descent {

void LocalExtreme::set_edge_nodes(math::Polynomial<double> const& fourth_degree_approximation,
    std::array<LocalExtreme::const_iterator, 7> const& samples, int i0, int i1, int i2)
{
  int left_index = 0;
  int right_index = 0;
  for (int j = 1; j <= 2; ++j)
  {
    if (samples[j]->w() < samples[left_index]->w())
      left_index = j;
    if (samples[j]->w() > samples[right_index]->w())
      right_index = j;
  }
  edge_nodes_[0] = samples[left_index];
  edge_nodes_[1] = samples[right_index];

  fourth_degree_approximation_ = fourth_degree_approximation;
}

#ifdef CWDEBUG
void LocalExtreme::debug_print_label(char const* left_or_right, const_iterator neighbor) const
{
  Dout(dc::notice, "Setting " << left_or_right << " neighbor of \"" << label() << "\" to \"" << neighbor->local_extreme().label() << "\".");
}
#endif

} // namespace gradient_descent
