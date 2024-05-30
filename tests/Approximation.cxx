#include "sys.h"
#include "Approximation.h"

namespace gradient_descent {

ScaleUpdate Approximation::add(Sample const* current, bool update_scale_only, bool current_is_replacement)
{
  DoutEntering(dc::notice, "Approximation::add(" << *current << ", " << std::boolalpha << update_scale_only << ", " <<
      current_is_replacement << ")");

  bool already_had_two_relevant_samples = number_of_relevant_samples_ == 2;
  if (current_is_replacement)
  {
    // The sample that was replaced in the history should already be in this approximation.
    ASSERT(number_of_relevant_samples_ > 0);
    relevant_samples_[current_index_] = current;
  }
  else
  {
    current_index_ = 1 - current_index_;
    relevant_samples_[current_index_] = current;
    number_of_relevant_samples_ = already_had_two_relevant_samples ? 2 : number_of_relevant_samples_ + 1;
  }

  if (!update_scale_only)
  {
    // LocalExtreme should only store parabolas of which the vertex was already determined;
    // hence changing it (by adding a new sample) doesn't make sense.
    ASSERT(!is_extreme_);

    if (number_of_relevant_samples_ == 1)
    {
      // If we have just one point, then the approximation is a linear function:
      //
      // A(w) = current->Lw() + L'(w) * (w - current->w())
      parabola_[2] = 0.0;
      parabola_[1] = current->dLdw();
      parabola_[0] = current->Lw() - current->dLdw() * current->w();
    }
    else
    {
      Sample const* prev = relevant_samples_[1 - current_index_];

      // The theory of this approach is described here:
      // https://math.stackexchange.com/a/4924305/489074
      double w0 = prev->w();
      double Lw0 = prev->Lw();
      double dLdw0 = prev->dLdw();
      double w1 = current->w();
      double Lw1 = current->Lw();
      double dLdw1 = current->dLdw();
      double delta_w = w0 - w1;
      double delta_Lw = Lw0 - Lw1;
      double delta_dLdw = dLdw0 - dLdw1;
      double sum_dLdw = dLdw0 + dLdw1;
      double product_Lw = dLdw0 * dLdw1;

      parabola_[2] = product_Lw < 0.0 ?
        (delta_Lw * sum_dLdw - 2.0 * delta_w * product_Lw) / (utils::square(delta_w) * delta_dLdw) :
        delta_Lw * delta_dLdw / (utils::square(delta_w) * sum_dLdw);

      parabola_[1] = delta_Lw / delta_w - (w0 + w1) * parabola_[2];
      parabola_[0] = (w0 * Lw1 - w1 * Lw0) / delta_w + w0 * w1 * parabola_[2];

      Dout(dc::notice, "Parabola fit through (" << w0 << ", " << Lw0 << ") with derivative " << dLdw0);
      Dout(dc::notice, "                 and (" << w1 << ", " << Lw1 << ") with derivative " << dLdw1);
      Dout(dc::notice, "    P(x) = " << parabola_ << " --> P(" << w0 << ") = " << parabola_(w0) << " and P(" << w1 << ") = " << parabola_(w1));
      ASSERT(utils::almost_equal(Lw0, parabola_(w0), 1e-6) && utils::almost_equal(Lw1, parabola_(w1), 1e-6));
      Dout(dc::notice, "    P'(" << w0 << ") = " << parabola_.derivative(w0) << " (requested: " << dLdw0 << ")");
      Dout(dc::notice, "    P'(" << w1 << ") = " << parabola_.derivative(w1) << " (requested: " << dLdw1 << ")");
    }
  }

  ScaleUpdate result = ScaleUpdate::first_sample;
  // Don't update the scale when current was just replaced. In that case we return first_sample,
  // but the caller passed current_is_replacement, so it should know to ignore that.
  if (!current_is_replacement && number_of_relevant_samples_ == 2)
    result = parabola_scale_.update(relevant_samples_, current_index_, parabola_, already_had_two_relevant_samples);
  return result;
}

} // namespace gradient_descent
