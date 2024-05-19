#include "sys.h"
#include "Approximation.h"

namespace gradient_descent {

std::string to_string(HorizontalDirection hdirection)
{
  switch (hdirection)
  {
    AI_CASE_RETURN(unknown_horizontal_direction);
    AI_CASE_RETURN(left);
    AI_CASE_RETURN(right);
  }
  AI_NEVER_REACHED
}

ScaleUpdate Approximation::add(Sample const* current, bool update_scale_only, bool current_is_replacement)
{
  Dout(dc::notice, "Approximation::add(" << *current << ", " << std::boolalpha << update_scale_only << ", " << current_is_replacement << ")");

  bool already_had_two_relevant_samples = number_of_relevant_samples_ == 2;
  if (current_is_replacement)
  {
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
    ASSERT(!is_part_of_extreme_);

    if (number_of_relevant_samples_ == 1)
    {
      // If we have just one point, then the approximation is a linear function:
      //
      // A(w) = coef[0] + L'(w) w
      parabola_[2] = 0.0;
      parabola_[1] = current->dLdw();
    }
    else
    {
      Sample const* prev = relevant_samples_[1 - current_index_];

      // If we have at least two points, the approximation is a parabola:
      //
      //   A(w) = approximation[0] + b w + c w²
      //
      // for which we determined the value of the derivative at two points:
      //
      //   L'(w₀) = b + 2c w₀
      //   L'(w₁) = b + 2c w₁
      //
      // In matrix form:
      //
      //   ⎡1 2w₀⎤ ⎡b⎤   ⎡L'(w₀)⎤
      //   ⎣1 2w₁⎦ ⎣c⎦ = ⎣L'(w₁)⎦
      //
      // from which follows
      //
      // ⎡b⎤        1     ⎡2w₁ -2w₀⎤⎡L'(w₀)⎤        1      ⎡2w₁ L'(w₀) - 2w₀ L'(w₁)⎤
      // ⎣c⎦ = ---------- ⎣-1   1  ⎦⎣L'(w₁)⎦ = ----------- ⎣    L'(w₁) -     L'(w₀)⎦
      //        2w₁ - 2w₀                      2 (w₁ - w₀)

      double inverse_det = 0.5 / (current->w() - prev->w());
      parabola_[2] = inverse_det * (current->dLdw() - prev->dLdw());
      // See https://math.stackexchange.com/questions/4913175
      parabola_[1] = inverse_det * (2.0 * current->Lw() - 2.0 * prev->Lw() + (prev->dLdw() - current->dLdw()) * (current->w() + prev->w()));
    }
    parabola_[0] = 0.0;
    parabola_[0] = current->Lw() - parabola_(current->w());
  }

  ScaleUpdate result = ScaleUpdate::first_sample;
  // Don't update the scale when current was just replaced. In that case we return first_sample,
  // but the caller passed current_is_replacement, so it should know to ignore that.
  if (!current_is_replacement && number_of_relevant_samples_ == 2)
    result = parabola_scale_.update(relevant_samples_, current_index_, parabola_, already_had_two_relevant_samples);
  return result;
}

} // namespace gradient_descent
