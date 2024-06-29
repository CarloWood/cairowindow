#include "sys.h"
#include "Approximation.h"

namespace gradient_descent {

ScaleUpdate Approximation::add(Sample const* current, bool current_is_replacement)
{
  DoutEntering(dc::notice|continued_cf, "Approximation::add(" << *current << ", " <<
      std::boolalpha << current_is_replacement << ") --> ");

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

  // Make sure that the samples are ordered, if we have two of them.
  if (number_of_relevant_samples_ == 2 && relevant_samples_[0]->w() > relevant_samples_[1]->w())
  {
    std::swap(relevant_samples_[0], relevant_samples_[1]);
    current_index_ = 1 - current_index_;
  }

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
    double delta_w_squared = utils::square(delta_w);
    double delta_Lw = Lw0 - Lw1;
    double delta_dLdw = dLdw0 - dLdw1;
    double sum_w = w0 + w1;
    double sum_dLdw = dLdw0 + dLdw1;
    double product_Lw = dLdw0 * dLdw1;

    parabola_[2] = product_Lw < 0.0 ?
      (delta_Lw * sum_dLdw - 2.0 * delta_w * product_Lw) / (delta_w_squared * delta_dLdw) :
      delta_Lw * delta_dLdw / (delta_w_squared * sum_dLdw);

    // If a parabola fit is not calculatable then treat it as a straight line through just the current sample.
    if (std::isnan(parabola_[2]) || std::isinf(parabola_[2]) ||
        std::abs(parabola_[2]) < Scale::epsilon || std::abs(delta_dLdw) < Scale::epsilon || std::abs(sum_dLdw) < Scale::epsilon)
    {
      parabola_[2] = 0.0;
      number_of_relevant_samples_ = 1;        // Ignore the oldest sample.
      parabola_scale_.reset();
      Dout(dc::finish, "ScaleUpdate::first_sample");
      return ScaleUpdate::first_sample;
    }

    parabola_[1] = delta_Lw / delta_w - sum_w * parabola_[2];
    parabola_[0] = (w0 * Lw1 - w1 * Lw0) / delta_w + w0 * w1 * parabola_[2];

    Dout(dc::notice, "Parabola fit through (" << w0 << ", " << Lw0 << ") with derivative " << dLdw0);
    Dout(dc::notice, "                 and (" << w1 << ", " << Lw1 << ") with derivative " << dLdw1);
    Dout(dc::notice, "    P(x) = " << parabola_ << " --> P(" << w0 << ") = " << parabola_(w0) << " and P(" << w1 << ") = " << parabola_(w1));
    Dout(dc::notice, "    P'(" << w0 << ") = " << parabola_.derivative(w0) << " (requested: " << dLdw0 << ")");
    Dout(dc::notice, "    P'(" << w1 << ") = " << parabola_.derivative(w1) << " (requested: " << dLdw1 << ")");

    // The theory of this approach is described here:
    // https://math.stackexchange.com/a/4926903/489074
    double d = (-2.0 * delta_Lw + delta_w * sum_dLdw) / (delta_w_squared * delta_w);
    double c = delta_dLdw / (2.0 * delta_w) - 1.5 * sum_w * d;
    double b = (w0 * dLdw1 - w1 * dLdw0) / delta_w + 3.0 * w0 * w1 * d;

    cubic_[3] = d;
    cubic_[2] = c;
    cubic_[1] = b;
    cubic_[0] = 0.0;
    cubic_[0] = Lw0 - cubic_(w0);

    Dout(dc::notice, "cubic = " << cubic_);
  }

  ScaleUpdate result = ScaleUpdate::first_sample;
  // Don't update the scale when current was just replaced. In that case we return first_sample,
  // but the caller passed current_is_replacement, so it should know to ignore that.
  if (!current_is_replacement && number_of_relevant_samples_ == 2)
    result = parabola_scale_.update(relevant_samples_, current_index_, parabola_, already_had_two_relevant_samples);

  Dout(dc::finish, result);
  return result;
}

ScaleUpdate Approximation::update_scale(Sample const& current)
{
  DoutEntering(dc::notice|continued_cf, "Approximation::update_scale(" << current << ") --> ");

  // This function should only be called on Approximation's that are part of a LocalExtreme.
  ASSERT(is_extreme_ && number_of_relevant_samples_ == 2);

  std::array<Sample const*, 2> samples = {{ &current, relevant_samples_[current_index_] }};
  ScaleUpdate result = parabola_scale_.update(samples, 0, parabola_, true);

  Dout(dc::finish, result);
  return result;
}

// find_extreme
//
// Return an extreme, if any, of the cubic approximation.
//
// Let `sl` be the left-most sample and `rl` the right-most sample,
// aka (relevant_samples_[0]) and relevant_samples_[1] respectively,
// because that is how relevant_samples_ is sorted).
//
// Input:
//  hdirection:
//         left: looking for an extreme left of `sr`.
//        right: looking for an extreme right of `sl`.
//    undecided: looking for an extreme, anywhere.
//  extreme:
//           up: looking for a maximum.
//         down: looking for a minimum.
//      unknown: looking for any extreme; the preference will depend on the value of hdirection.
//               If hdirection is left, the extreme nearest to, but left of, sr is returned.
//               If hdirection is right, the extreme nearest to, but right of, sl is returned.
//               if hdirection is undecided too, the prefered extreme is the one nearest to the center.
//
// Output:
//  hdirection:
//         left: the requested extreme is left of `left`.
//        right: the requested extreme is right of `right`.
//    inbetween: the requested extreme is in between `left` and `right`.
//  extreme:
//           up: the extreme is a maximum.
//         down: the extreme is a minimum.
//      unknown: there is no extreme (of the requested vdirection type),
//               in this case, if hdirection was undecided, it is set to
//               point in the direction where we go downhill. Otherwise
//               it is unchanged.
Weight Approximation::find_extreme(HorizontalDirection& hdirection, VerticalDirection& extreme) const
{
  DoutEntering(dc::notice, "Approximation::find_extreme(" << hdirection << ", " << extreme << ")");
  // Otherwise the cubic approximation isn't even valid.
  ASSERT(number_of_relevant_samples_ == 2);

  double a = cubic_[0];
  double b = cubic_[1];
  double c = cubic_[2];
  double d = cubic_[3];

  double D = utils::square(c) - 3.0 * b * d;

  // If D is equal to zero we have a point where the derivative is zero, but that isn't a maximum or minimum.
  if (D <= 0.0)
  {
    extreme = VerticalDirection::unknown;
    if (hdirection == HorizontalDirection::undecided)
      hdirection = (d > 0.0 || (d == 0.0 && b > 0.0)) ? HorizontalDirection::left : HorizontalDirection::right;
    return {};
  }

  double sqrt_D = std::sqrt(D);
  double minimum = (-c + sqrt_D) / (3.0 * d);
  double maximum = (-c - sqrt_D) / (3.0 * d);

  double wl = relevant_samples_[0]->w();
  double wr = relevant_samples_[1]->w();

  if (extreme == VerticalDirection::unknown)
  {
    //      unknown: looking for any extreme; the preference will depend on the value of hdirection.
    //               If hdirection is left, the extreme nearest to, but left of, sr is returned.
    //               If hdirection is right, the extreme nearest to, but right of, sl is returned.
    //               if hdirection is undecided too, the prefered extreme is the one nearest to the center.
    if (hdirection == HorizontalDirection::undecided)
    {
      double center = 0.5 * (wl + wr);
      extreme = (std::abs(minimum - center) > std::abs(maximum - center)) ? VerticalDirection::up : VerticalDirection::down;
    }
    else
    {
      ASSERT(hdirection == HorizontalDirection::left || hdirection == HorizontalDirection::right);
      double reference = (hdirection == HorizontalDirection::left) ? wr : wl;
      double rel_pos_min = static_cast<int>(hdirection) * (minimum - reference);
      double rel_pos_max = static_cast<int>(hdirection) * (maximum - reference);
      if (rel_pos_min < 0.0 && rel_pos_max < 0.0)
      {
        // Both extremes are on the wrong side of reference.
        extreme = VerticalDirection::unknown;
        return {};
      }
      if (rel_pos_min < 0.0)
        extreme = VerticalDirection::up;        // The minimum is on the wrong side of reference.
      else if (rel_pos_max < 0.0)
        extreme = VerticalDirection::down;      // The maximum is on the wrong side of reference.
      else
        extreme = (rel_pos_max < rel_pos_min) ? VerticalDirection::up : VerticalDirection::down;
    }
  }

  double result = (extreme == VerticalDirection::up) ? maximum : minimum;

  //         left: looking for an extreme left of `sr`.
  //        right: looking for an extreme right of `sl`.
  if ((hdirection == HorizontalDirection::left && result >= wr) ||
      (hdirection == HorizontalDirection::right && result <= wl))
  {
    extreme = VerticalDirection::unknown;
    return {};
  }

  //         left: the requested extreme is left of `left`.
  //        right: the requested extreme is right of `right`.
  //    inbetween: the requested extreme is in between `left` and `right`.
  if (result < wl)
    hdirection = HorizontalDirection::left;
  else if (result > wr)
    hdirection = HorizontalDirection::right;
  else
    hdirection = HorizontalDirection::inbetween;

  return result;
}

} // namespace gradient_descent
